#include "fvCFD.H"
#include "polyMesh.H"
#include "IOstreams.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "Pstream.H"  // parallel communication
#include "argList.H"
#include <cstdlib>

int main(int argc, char *argv[])
{
    // -- Select probe point option
    argList::addOption
    (
        "point",
        "vector",
        "create cutting planes normal to the given direction <vector> - eg, '(1 0 0)'"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // -- go through all times
    instantList timeDirs = timeSelector::select0(runTime, args);
    
    // -- read probe vector
    vector X_physVec(0,0,1);
    if (args.optionReadIfPresent("point", X_physVec))
    {
        Info<< "Point for probe: " << X_physVec << endl;
    }
    else
    {
        Info<< "Default point for probe: " << X_physVec << endl;
    }
    point X_phys = point(X_physVec);

    for (int i = 0; i < timeDirs.size(); i++)
    {
        runTime.setTime(timeDirs[i], timeDirs.size()-1);   

        // -- Parameters for Newton
        int maxIter = 30;
        scalar tol = 1e-9;

        // -----------------------------
        // Read fields
        // -----------------------------
        volVectorField D
        (
            IOobject(
                "D", 
                runTime.timeName(), 
                mesh,
                IOobject::MUST_READ, 
                IOobject::NO_WRITE
            ), 
            mesh
        );
        volTensorField gradD
        (
            IOobject
            (
                "grad(D)", 
                runTime.timeName(), 
                mesh,
                IOobject::MUST_READ, 
                IOobject::NO_WRITE
            ), 
            mesh
        );

        volScalarField T
        (
            IOobject(
                "T", 
                runTime.timeName(), 
                mesh,
                IOobject::MUST_READ, 
                IOobject::NO_WRITE
            ), 
            mesh
        );

        // -----------------------------
        // Newton iteration to find X_ref
        // -----------------------------
        // bool probeFound = false;
        bool probeFoundThisProc = false;
        point X_ref = X_phys;

        for (int iter=0; iter<maxIter; iter++)
        {
            label cellI = mesh.findCell(X_ref);
            probeFoundThisProc = (cellI >= 0);

            vector Di(0,0,0);
            tensor gradDi(tensor::zero);
            autoPtr<interpolation<vector>> Dinterp(nullptr);
            autoPtr<interpolation<tensor>> gradDinterp(nullptr);

            if (mesh.nCells() > 0)  // construct on all ranks
            {
                dictionary interpolationDict = mesh.solutionDict().subDict("interpolationSchemes");
                Dinterp = interpolation<vector>::New(interpolationDict, D);
                gradDinterp = interpolation<tensor>::New(interpolationDict, gradD);
            }
            if (probeFoundThisProc)  // this processor owns the cell
            {
                Di = Dinterp->interpolate(X_ref, cellI);
                gradDi = gradDinterp->interpolate(X_ref, cellI);
            }

            tensor J = tensor::I + gradDi;
            vector F = X_ref + Di - X_phys;
            vector deltaX = inv(J) & F;

            if (!probeFoundThisProc) deltaX = vector(-1,-1,-1);

            vector globalDeltaX = deltaX;
            Foam::reduce(globalDeltaX, Foam::maxOp<vector>());
            if (mag(globalDeltaX) < tol)
            {
                break;
            }
            X_ref -= globalDeltaX;
        }

        label cellI = mesh.findCell(X_ref);
        probeFoundThisProc = (cellI >= 0);
        autoPtr<interpolation<scalar>> Tinterp(nullptr);
        scalar Ti = -1;
        if (mesh.nCells() > 0)  // construct on all ranks
        {
            dictionary interpolationDict = mesh.solutionDict().subDict("interpolationSchemes");
            Tinterp = interpolation<scalar>::New(interpolationDict, T);
        }
        if (probeFoundThisProc)  // this processor owns the cell
        {
            Ti = Tinterp->interpolate(X_ref, cellI);
        }
        scalar globalT = Ti;
        Foam::reduce(globalT, Foam::maxOp<scalar>());
        Info << "Time = " << runTime.timeName() << "; T = " << globalT << endl;
    }
    Info << "End" <<endl;
    return 0;
}