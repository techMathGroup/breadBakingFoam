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

        // -- CO2 bookkeeping (initialized on first timestep)
        static bool co2Init = false;
        static double initialBreadMass = 0.0;  // total mass (solid + liquid + gas) at t0
        static double cumCO2 = 0.0;            // cumulative CO2 produced [kg]

        physics().setDeltaT(runTime);

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        volScalarField invT = 1/T;

        const double dtVal = runTime.deltaT().value();

        if (!co2Init)
        {
            const dimensionedScalar breadMass =
                fvc::domainIntegrate
                (
                    J *
                    (
                        alphaS * rhoS
                      + alphaL * rhoL
                      + (scalar(1) - alphaS - alphaL) * rhoG
                    )
                );

            initialBreadMass = breadMass.value();
            co2Init = true;
        }

        while (pimple.loop())
        {
            // -- update moisture content for evaporation source calculation
            moisture = (alphaL * rhoL) / (alphaS * rhoS);

            // -- evaporation source calculation 
            #include "compEvS.H"

            // -- update of the effective diffusivity, heat conductivity and permeability
            #include "compEffProps.H"

            // -- fermentation source
            mCO2 = R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2)) * rhoS / dummyRho * alphaS;
            
            // -- molar amount change
            nCO2 = fvc::ddt(alphaG, rhoG, J)/Mg;
            nCO2.correctBoundaryConditions();

            // -- calculate pre-coefficients for flux calculations
            jGTilda = - rhoG * permGLViscG * Finv.T();
            jDVTilda = - rhoG * DEff * Finv.T();

            // -- conservation of the water vapor
            #include "concEqC5.H"
            #include "concEqV5.H"
            omegaAir = 1 - omegaV - omegaC;
            omegaAir.correctBoundaryConditions();
            // #include "concEqAir.H"

            // -- gas molar mass
            Mg = 1 / (omegaV / molMV + omegaC / molMC + (1 - omegaV - omegaC) / molMAir);
            Mg.correctBoundaryConditions();

            // -- energy conservation
            #include "EEqn2.H"

            // -- crust formation
            #include "crustEq.H"
            
            // -- if deformation allowed
            if (withDeformation == 1)
            {
                for (int i = 0; i < 1; ++i)
                {
                    // -- dilatation laws --> solids4foam with custom visco-elastic model
                    physics->evolve();
                    
                    // -- deformation gradient
                    F = I + physics->gradD().T();

                    // -- deformation gradient inverse
                    Finv = inv(F);
                    Finv.correctBoundaryConditions();

                    // -- Jacobian of the deformation gradient
                    J = det(F);
                    J.correctBoundaryConditions();
                }
            }

            
            // -- solid mass conservation
            #include "alphaSEq.H"
            
            // -- liquid water conservation
            #include "phiLEq.H"

            alphaG = 1 - alphaS - alphaL;
            alphaG.correctBoundaryConditions();

            #include "concEqG6.H"

            rhoG = Mg / univR / T * pG;
            rhoG.correctBoundaryConditions();

            // -- basic log
            if (debug >= 1)
            {
                Info << "phiL   : res: " << phiLResidual << " Min (phiL): " << min(alphaL).value() << ", max (phiL): " << max(alphaL).value() << "." << endl;
                Info << "pG     : res: " << pResidual    << " Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
                Info << "T      : res: " << TResidual    << " Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
                Info << "omV    : res: " << omegaVResidual << " Min (omegaV): " << min(omegaV).value() << ", max (omegaV): " << max(omegaV).value() << "." << endl;
                Info << "omC    : res: " << omegaCResidual << " Min (omegaC): " << min(omegaC).value() << ", max (omegaC): " << max(omegaC).value() << "." << endl;
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
