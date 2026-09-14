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
        Info << "Time = " << runTime.timeName() << endl;
        runTime.setTime(timeDirs[i], timeDirs.size()-1);
        volVectorField D
        (
            IOobject
                (
                // "moisturePostProcess", 
                "D", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField rhoG
        (
            IOobject
                (
                "rhoG", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // Create output directory and file once per time step
        fileName outputDir = runTime.path();
        mkDir(outputDir/"../Shape");
        fileName outputFile = outputDir/"../Shape/D_values_" + runTime.timeName() + ".dat";
        bool intro = false;

        // Initialize file (clear it on first write)
        if (Pstream::master())
        {
            std::ofstream initFile(outputFile.c_str(), std::ios::trunc);
            initFile.close();
        }

        // surfaceVectorField Cf = mesh.Cf();

        forAll(D.boundaryField(), patchI)
        {
            if (rhoG.boundaryField()[patchI].type() == "calculated")
            {
                vectorField& DP = D.boundaryFieldRef()[patchI];
                const vectorField& Cf = mesh.boundary()[patchI].Cf();
                DP = DP + Cf;
                
                labelList sizes(Pstream::nProcs());
                sizes[Pstream::myProcNo()] = DP.size();
                Pstream::gatherList(sizes);
                Pstream::scatterList(sizes);
                
                if (Pstream::master())
                {
                    Info << "Writing patch " << patchI << " to: " << outputFile << endl;
                    
                    std::ofstream os(outputFile.c_str(), std::ios::app);

                    if (!intro)
                    {
                        os << "x" << "\t" << "y" << "\t" << "z" << "\n";
                        intro = true;
                    }
                    // os << "x" << "\t" << "y" << "\t" << "z" << "\n";
                    
                    if (!os.good())
                    {
                        FatalErrorInFunction
                            << "Cannot open file: " << outputFile
                            << exit(FatalError);
                    }
                    
                    forAll(sizes, procI)
                    {
                        vectorField procDP = DP;
                        if (procI != Pstream::myProcNo())
                        {
                            procDP.resize(sizes[procI]);
                            IPstream fromProc(UPstream::commsTypes::scheduled, procI);
                            fromProc >> procDP;
                        }
                        
                        forAll(procDP, faceI)
                        {
                            os << procDP[faceI][0] << "\t" << procDP[faceI][1] << "\t" << procDP[faceI][2] << "\n";
                        }
                    }
                    os.close();
                }
                else
                {
                    OPstream toMaster(UPstream::commsTypes::scheduled, Pstream::masterNo());
                    toMaster << DP;
                }
            }
        }
        rhoG.correctBoundaryConditions();
        // volScalarField J
        // (
        //     IOobject
        //         (
        //         "J", 
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::MUST_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh
        // );

        // Compute local partial sums
        // scalar localSum = gSum(moisture.internalField() * J.internalField() * mesh.V().field());
        // scalar localSum = sum(moisture.internalField()  * J.internalField() * mesh.V().field());
        // scalar localSum = sum(moisture.internalField()   * mesh.V().field());


        // scalar totalVol  = sum(mesh.V().field() );

        // Parallel reduction
        // scalar globalSum = localSum;
        // Foam::reduce(globalSum, Foam::sumOp<scalar>());
        // scalar globalVol  = totalVol;
        // Foam::reduce(globalVol, Foam::sumOp<scalar>());

        // scalar avg = globalSum / globalVol;
        // scalar avgB = globalSumB / globalVol;

        // Info << "Time = " << runTime.timeName() << "; Moisture average = " << avg << endl;

    }
    Info << "End" << endl;
    return 0;
}


// ************************************************************************* //
