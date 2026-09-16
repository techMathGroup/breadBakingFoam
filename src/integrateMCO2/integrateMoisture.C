/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    caclPressDerOnSlices

Description
    Calculates and prints average pressures over slices along some given
    coordinate. In user specified box
     

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "sampledPlane.H"
// #include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // #include "setRootCase.H"
    // #include "createTime.H"
    // instantList timeDirs = timeSelector::select0(runTime, args);
    // #include "createNamedMesh.H"

    // scalar V = 0;
    // // Info <<timeDirs.size()<<endl;
    // for (int i = 0; i < timeDirs.size(); i++)
    // {
    //     runTime.setTime(timeDirs[i], timeDirs.size()-1);
    //     // Info<< "Time = " << runTime.timeName() << endl;
    //     mesh.readUpdate();

    //     volScalarField moisture
    //     (
    //         IOobject
    //             (
    //             "moisture", 
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::MUST_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh
    //     );

    //     volScalarField J
    //     (
    //         IOobject
    //             (
    //             "J", 
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::MUST_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh
    //     );

    //     // -- integration (sum c_s*V in cellZone)
    //     scalar integral(0);
    //     scalar VTu(0);

    //     // // -- loop over all cells in mesh    
    //     forAll(mesh.cells(), celli)
    //     {
    //         integral += mesh.V()[celli] * moisture[celli] * J[celli];
    //         VTu += mesh.V()[celli] * J[celli];
    //     }
    //     V = VTu;
    //     Pout << "Time = " << runTime.timeName() << "; Moisture average = " << integral/VTu << endl;
    // }

    // Pout<< "End" << V << endl;

    // return 0;

    // argList::validArgs.append("case");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    for (int i = 0; i < timeDirs.size(); i++)
    {
        runTime.setTime(timeDirs[i], timeDirs.size()-1);
        volScalarField mCO2
        (
            IOobject
                (
                // "moisturePostProcess", 
                "mCO2", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField J
        (
            IOobject
                (
                "J", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField mCO2Dis
        (
            IOobject
                (
                "mCO2Dis", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField evRate
        (
            IOobject
                (
                "evRate", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // volScalarField alphaG
        // (
        //     IOobject
        //         (
        //         "alphaG", 
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::MUST_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh
        // );

        // volScalarField rhoG
        // (
        //     IOobject
        //         (
        //         "rhoG", 
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::MUST_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh
        // );

        // volScalarField omegaV
        // (
        //     IOobject
        //         (
        //         "omegaV", 
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::MUST_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh
        // );

        // Compute local partial sums
        scalar localSumMCO2 =    sum(J.internalField() * mCO2.internalField() * mesh.V().field());
        scalar localSumMCO2Dis = sum(J.internalField() * mCO2Dis.internalField() * mesh.V().field());
        scalar localSumEvRate =  sum(J.internalField() * evRate.internalField() * mesh.V().field());


        scalar totalVol  = sum(mesh.V().field() * J.internalField());

        // Parallel reduction
        // scalar globalSum = returnReduce(localSum, sumOp<scalar>());
        // scalar globalVol  = returnReduce(totalVol, sumOp<scalar>());
        scalar globalSumMCO2 = localSumMCO2;
        Foam::reduce(globalSumMCO2, Foam::sumOp<scalar>());
        scalar globalSumMCO2Dis = localSumMCO2Dis;
        Foam::reduce(globalSumMCO2Dis, Foam::sumOp<scalar>());
        scalar globalSumEvRate = localSumEvRate;
        Foam::reduce(globalSumEvRate, Foam::sumOp<scalar>());
        scalar globalVol  = totalVol;
        Foam::reduce(globalVol, Foam::sumOp<scalar>());
        // rhoD.correctBoundaryConditions();

        // scalar avg = globalSum / globalVol;
        // scalar avg = globalSum;
        // scalar avgB = globalSumB / globalVol;

        // Info << "Time = " << runTime.timeName() << "; rhoD average = " << avg << endl;
        // Info << "Time = " << runTime.timeName() << "; weight = " << globalSum << "total volume = " << globalVol << endl;
        // Info << "Time = " << runTime.timeName() << "; " << globalSumMCO2 / totalVol << "; " << globalSumMCO2Dis / totalVol << "; " << globalSumEvRate / totalVol << endl;
        Info << "Time = " << runTime.timeName() << " " << globalSumMCO2  << " " << globalSumMCO2Dis  << " " << globalSumEvRate << " " << globalVol << endl;

    }
    Info << "End" << endl;
    return 0;
}


// ************************************************************************* //
