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

        physics().setDeltaT(runTime);

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        volScalarField invT = 1/T;

        dimensionedScalar dTime = runTime.time().deltaT();

        while (pimple.loop())
        {
            // -- if deformation allowed
            if (withDeformation == 1)
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

            // -- solid mass conservation
            #include "alphaSEq.H"
            
            // -- liquid water conservation
            #include "phiLEq.H"

            // -- update of the effective diffusivity, heat conductivity and permeability
            #include "compEffProps.H"
            
            // -- gas conservation 
            #include "concEqG6.H"
            
            // -- gas density update
            rhoG = pG * Mg / univR / T;
            rhoG.correctBoundaryConditions();

            // -- conservation of the water vapor
            #include "concEqV5.H"
            
            // -- energy conservation
            #include "EEqn2.H"

            // -- CO2 source by fermentation (Zhang Transport Processes and Large Deformation During Baking of Bread 2005)
            mCO2 = R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2)) * rhoS / dimensionedScalar("dummyRho", dimMass/dimVolume, 1) * alphaS;

            // -- moisture content -- kgLwater / kgS
            moisture = (alphaL * rhoL + (1 - alphaL - alphaS) * rhoG * omegaV) / (alphaS * rhoS);

            // -- evaporation source calculation 
            #include "compEvS.H"

            // -- basic log
            if (debug >= 1)
            {
                Info << "phiL   : res: " << phiLResidual << " Min (phiL): " << min(alphaL).value() << ", max (phiL): " << max(alphaL).value() << "." << endl;
                Info << "pG     : res: " << pResidual    << " Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
                Info << "T      : res: " << TResidual    << " Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
                Info << "omV    : res: " << omegaVResidual << " Min (omegaV): " << min(omegaV).value() << ", max (omegaV): " << max(omegaV).value() << "." << endl;
                Info << "Min (J): " << min(J).value() << ", max (J): " << max(J).value() << "." << endl;
                Info << "Min (permGLViscG): " << min(permGLViscG).value() << ", max (permGLViscG): " << max(permGLViscG) << "." << endl;
                Info << endl;
            }
        }

        physics().updateTotalFields();

        runTime.printExecutionTime(Info);
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
