/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    breadBakingFoam

Group
    --

Description
    Transient solver for bread baking

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
// #include "fvOptions.H"
#include "physicsModel.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for compressible turbulent flow,"
        "with implicit or explicit porosity treatment and optional sources."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    autoPtr<physicsModel> physics = physicsModel::New(runTime);

    #include "readTransportProperties.H"
    #include "createControlMy.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run()){

        SolverPerformance<scalar>::debug = 0;
        SolverPerformance<vector>::debug = 0;

        physics().setDeltaT(runTime);

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        pG.storePrevIter();
        J.storePrevIter();
        omegaV.storePrevIter();
        omegaC.storePrevIter();
        omegaAir.storePrevIter();
        Finv.storePrevIter();
        jDVTilda.storePrevIter();
        jDCTilda.storePrevIter();
        jDATilda.storePrevIter();
        int iter = 0;

        while (pimple.loop())
        {
            iter++;
            if (withDeformation == 1)
            {
                // -- deformation laws --> solids4foam with custom visco-elastic model
                physics->evolve();
                
                // -- deformation gradient
                F = physics->mesh().lookupObject<volTensorField>("F");
                F.correctBoundaryConditions();

                // -- deformation gradient inverse
                Finv = physics->mesh().lookupObject<volTensorField>("Finv");
                Finv.correctBoundaryConditions();

                // -- Jacobian of the deformation gradient
                J = physics->mesh().lookupObject<volScalarField>("J");
                J.correctBoundaryConditions();
            }
            
            // -- compute volumetric fractions
            alphaD = alphaD0 / J;
            alphaD.correctBoundaryConditions();

            alphaG = 1.0 - alphaD;
            alphaG.correctBoundaryConditions();

            int nIter = 1;
            if (pimple.finalIter())
            {
                nIter = nItersAfter;
            }

            for (int i = 0; i < nIter; ++i)
            {

                // -- evaporation source calculation 
                #include "compEvRate.H"

                // -- update of the effective diffusivity, heat conductivity and permeability
                #include "compEffProps.H"
                
                // -- liquid water conservation
                #include "phiLEq.H"

                // -- dissolved CO2 conservation
                #include "disCO2Eqn.H"

                // -- fermentation source (kg/m3/s)
                mCO2 = R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2));
                mCO2.correctBoundaryConditions();
                
                // -- calculation of pre-coefficients for flux calculations
                jGTilda = - rhoG * gasDarcyMobility * Finv.T();
                jGTilda.correctBoundaryConditions();

                jDVTilda = - rhoG * DEffvM * Finv.T();
                jDVTilda.correctBoundaryConditions();
                jDCTilda = - rhoG * DEffcM * Finv.T();
                jDCTilda.correctBoundaryConditions();
                jDATilda = - rhoG * DEffaM * Finv.T();
                jDATilda.correctBoundaryConditions();
            
                // -- correction of diffusive flux
                jC = - ((jDVTilda & fvc::grad(omegaV)) + (jDCTilda & fvc::grad(omegaC)) + (jDATilda & fvc::grad(omegaAir)));
                jC.correctBoundaryConditions();

                jVE = omegaV * ((jGTilda & fvc::grad(pG)) + jC) + (jDVTilda & fvc::grad(omegaV)); 
                jCE = omegaC * ((jGTilda & fvc::grad(pG)) + jC) + (jDCTilda & fvc::grad(omegaC)); 
                jAE = omegaAir * ((jGTilda & fvc::grad(pG)) + jC) + (jDATilda & fvc::grad(omegaAir)); 
                jVE.correctBoundaryConditions();
                jCE.correctBoundaryConditions();
                jAE.correctBoundaryConditions();

                if (i % 2 == 0)
                {
                    // -- overall gas-phase balance
                    #include "concEqG.H"
                        
                    #include "EEqn.H"

                    // -- gas density calculation
                    rhoG = Mg / univR / T * pG;
                    rhoG.correctBoundaryConditions();
                    
                    // -- species equations
                    
                    #include "concEqC.H"
                    #include "concEqV.H"


                    // -- last species
                    omegaAir = 1.0 - omegaV - omegaC;
                    omegaAir.correctBoundaryConditions();

                }

                else
                {

                    // -- species equations
                    #include "concEqV.H"
                    #include "concEqC.H"

                    
                    // -- last species
                    omegaAir = 1.0 - omegaV - omegaC;
                    omegaAir.correctBoundaryConditions();

                    
                    // -- overall gas-phase balance
                    #include "concEqG.H"
                        
                    #include "EEqn.H"

                    // -- gas density calculation
                    rhoG = Mg / univR / T * pG;
                    rhoG.correctBoundaryConditions();
                }

                // -- last species
                omegaAir = 1.0 - omegaV - omegaC;
                omegaAir.correctBoundaryConditions();


                // -- gas properties calculation
                Mg = 1.0 / (omegaV / molMV + omegaC / molMC + omegaAir / molMAir);
                Mg.correctBoundaryConditions();

                // -- molar fractions
                yV = omegaV / molMV * Mg;
                yV.correctBoundaryConditions();
                yC = omegaC / molMC * Mg;
                yC.correctBoundaryConditions();
                yA = omegaAir / molMAir * Mg;
                yA.correctBoundaryConditions();

                // -- gas density calculation
                rhoG = Mg / univR / T * pG;
                rhoG.correctBoundaryConditions();

                if (pNIter == 0)
                {
                    break;
                }

            }

            

            // -- basic log
            if (debug >= 1)
            {
                Info << "phiL   : res: " << phiLResidual << " Min (moisture): " << min(moisture).value() << ", max (moisture): " << max(moisture).value() << "." << endl;
                Info << "pG     : res: " << pResidual    << " Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
                Info << "T      : res: " << TResidual    << " Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
                Info << "omV    : res: " << omegaVResidual << " Min (omegaV): " << min(omegaV).value() << ", max (omegaV): " << max(omegaV).value() << "." << endl;
                Info << "omC    : res: " << omegaCResidual << " Min (omegaC): " << min(omegaC).value() << ", max (omegaC): " << max(omegaC).value() << "." << endl;
                Info << "Min (alphaG): " << min(alphaG).value() << ", max (alphaG): " << max(alphaG).value() << "." << endl;
                Info << "Min (J): " << min(J).value() << ", max (J): " << max(J).value() << "." << endl;
                Info << "Min (gasDarcyMobility): " << min(gasDarcyMobility).value() << ", max (gasDarcyMobility): " << max(gasDarcyMobility) << "." << endl;
                Info << endl;
            }
        }

        physics().updateFields();

        physics().updateTotalFields();

        runTime.printExecutionTime(Info);
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
