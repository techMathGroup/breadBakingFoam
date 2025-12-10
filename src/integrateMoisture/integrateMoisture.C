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
        volScalarField moisture
        (
            IOobject
                (
                "moisture", 
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

        // Compute local partial sums
        scalar localSum = gSum(moisture.internalField() * J.internalField() * mesh.V().field());


        scalar totalVol  = gSum(mesh.V().field() * J.internalField());

        // Parallel reduction
        scalar globalSum = returnReduce(localSum, sumOp<scalar>());
        scalar globalVol  = returnReduce(totalVol, sumOp<scalar>());

        scalar avg = globalSum / globalVol;
        // scalar avgB = globalSumB / globalVol;

        Info << "Time = " << runTime.timeName() << "; Moisture average = " << avg << endl;

    }
    Info << "End" << endl;
    return 0;
}


// ************************************************************************* //
