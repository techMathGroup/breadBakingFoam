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
#include "fvOptions.H"
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

        #include "initCO2Calc.H"

        physics().setDeltaT(runTime);

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        pG.storePrevIter();
        jLCorr.storePrevIter();
        jLCorrIntoLiquid.storePrevIter();
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
                for (int i = 0; i < 1; ++i)
                {
                    // -- dilatation laws --> solids4foam with custom visco-elastic model
                    physics->evolve();
                    
                    // -- deformation gradient
                    F = I + physics->gradD().T();

                    // -- deformation gradient inverse
                    Finv = physics->mesh().lookupObject<volTensorField>("Finv");;
                    Finv.correctBoundaryConditions();

                    // -- Jacobian of the deformation gradient
                    J = physics->mesh().lookupObject<volScalarField>("J");
                    J.correctBoundaryConditions();
                }
            }

            alphaD = alphaD0 / J;
            alphaD.correctBoundaryConditions();

            alphaG = 1.0 - alphaD;
            alphaG.correctBoundaryConditions();

            // forAll(alphaD.boundaryField(), patchI)
            // {
            //     scalarField& alphaDP = alphaD.boundaryFieldRef()[patchI];
            //     scalarField& alphaGP = alphaG.boundaryFieldRef()[patchI];
            //     scalarField& JP = J.boundaryFieldRef()[patchI];

            //     forAll(alphaDP, faceI)
            //     {
            //         alphaDP[faceI] = alphaD0 / JP[faceI];
            //         alphaGP[faceI] = 1.0 - alphaDP[faceI];
            //     }
            // }



            scalar nConc = 1;
            if (pimple.finalIter())
            {
                nConc = 1;
            }

            for (int i = 0; i < nConc; i++)
            {

                // -- evaporation source calculation 
                #include "compEvS.H"

                // -- update of the effective diffusivity, heat conductivity and permeability
                #include "compEffProps.H"

                // -- solid mass conservation
                #include "alphaSEq.H"

                // -- closed-cell water vapor flux
                jL = - dKoeffLucas * ((Finv.T() & fvc::grad(awPSat)));
                // jL = - dKoeffLucas * ((Finv.T() & fvc::grad(pV)));
                jL.correctBoundaryConditions();
                
                // -- liquid water conservation
                #include "phiLEq.H"
                omegaL = 1.0 - omegaS;
                omegaL.correctBoundaryConditions();

                // -- update moisture content
                moisture = omegaL / omegaS;
                moisture.correctBoundaryConditions();

                // -- fermentation source
                // mCO2 = R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2)) * rhoS / dummyRho * alphaS;
                mCO2 = R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2)) * alphaD0 * rhoD * omegaS / dummyRho; 
                mCO2.correctBoundaryConditions();
                
                // -- calculation of pre-coefficients for flux calculations
                jGTilda = - rhoG * permGLViscG * Finv.T();
                jGTilda.correctBoundaryConditions();

                jDVTilda = - rhoG * DEffvM * Finv.T();
                jDVTilda.correctBoundaryConditions();
                jDCTilda = - rhoG * DEffcM * Finv.T();
                jDCTilda.correctBoundaryConditions();
                jDATilda = - rhoG * DEffaM * Finv.T();
                jDATilda.correctBoundaryConditions();
            
                // -- correction of diffusive flux
                jC = - 0*((jDVTilda & fvc::grad(omegaV)) + (jDCTilda & fvc::grad(omegaC)) + (jDATilda & fvc::grad(omegaAir)));
                jC.correctBoundaryConditions();

                jVE = omegaV * ((jGTilda & fvc::grad(pG)) + jC) + (jDVTilda & fvc::grad(omegaV)); 
                jCE = omegaC * ((jGTilda & fvc::grad(pG)) + jC) + (jDCTilda & fvc::grad(omegaC)); 
                // jAE = omegaAir * ((jGTilda & fvc::grad(pG)) + jC) + (jDATilda & fvc::grad(omegaAir)); 
                jVE.correctBoundaryConditions();
                jCE.correctBoundaryConditions();
                // jAE.correctBoundaryConditions();


                // -- overall gas-phase balance
                #include "concEqG6.H"

                // -- species equations
                #include "concEqV5.H"
                // #include "concEqC5.H"
                // #include "concEqA5.H"

                // -- last species
                // omegaAir = 1.0 - omegaV - omegaC;
                // omegaAir.correctBoundaryConditions();

                omegaC = 1.0 - omegaV;
                omegaC.correctBoundaryConditions();
                   
                #include "EEqn4.H"

                // -- gas density calculation
                rhoG = Mg / univR / T * pG;
                rhoG.correctBoundaryConditions();

                // -- gas properties calculation
                // Mg = 1.0 / (omegaV / molMV + omegaC / molMC + omegaAir / molMAir);
                Mg = 1.0 / (omegaV / molMV + omegaC / molMC);
                Mg.correctBoundaryConditions();

                // -- molar fractions
                yV = omegaV / molMV * Mg;
                yV.correctBoundaryConditions();
                yC = omegaC / molMC * Mg;
                yC.correctBoundaryConditions();
                // yA = omegaAir / molMAir * Mg;
                // yA.correctBoundaryConditions();

                // -- gas density calculation
                rhoG = Mg / univR / T * pG;
                rhoG.correctBoundaryConditions();

            }

            // -- basic log
            if (debug >= 1)
            {
                Info << "phiL   : res: " << phiLResidual << " Min (rhoD): " << min(rhoD).value() << ", max (rhoD): " << max(rhoD).value() << "." << endl;
                Info << "pG     : res: " << pResidual    << " Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
                Info << "T      : res: " << TResidual    << " Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
                Info << "omV    : res: " << omegaVResidual << " Min (omegaV): " << min(omegaV).value() << ", max (omegaV): " << max(omegaV).value() << "." << endl;
                Info << "omC    : res: " << omegaCResidual << " Min (omegaC): " << min(omegaC).value() << ", max (omegaC): " << max(omegaC).value() << "." << endl;
                Info << "Min (alphaG): " << min(alphaG).value() << ", max (alphaG): " << max(alphaG).value() << "." << endl;
                Info << "Min (J): " << min(J).value() << ", max (J): " << max(J).value() << "." << endl;
                Info << "Min (permGLViscG): " << min(permGLViscG).value() << ", max (permGLViscG): " << max(permGLViscG) << "." << endl;
                Info << endl;
            }
        }

        physics().updateFields();

        physics().updateTotalFields();

        // -- CO2 accumulation and optional log
        const dimensionedScalar mDotCO2 = fvc::domainIntegrate(J * mCO2);
        const double mDotCO2Val = mDotCO2.value();

        cumCO2 += mDotCO2Val * dtVal;

        if (debug >= 1)
        {
            Info<< "CO2 prod [kg/s]: " << mDotCO2Val
                << ", cum CO2 [kg]: " << cumCO2;

            if (initialBreadMass > SMALL)
            {
                Info<< ", CO2/bread ratio: " << cumCO2 / initialBreadMass;
            }

            Info<< endl;
        }

        runTime.printExecutionTime(Info);
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
