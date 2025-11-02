/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "viscoBread.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscoBread, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, viscoBread, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar viscoBread::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label viscoBread::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar viscoBread::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar viscoBread::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::viscoBread::makeJ()
{
    if (JPtr_.valid())
    {
        FatalErrorIn("void Foam::viscoBread::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JPtr_().oldTime();
}


Foam::volScalarField& Foam::viscoBread::J()
{
    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


void Foam::viscoBread::makeJf()
{
    if (JfPtr_.valid())
    {
        FatalErrorIn("void Foam::viscoBread::makeJf()")
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "lawJf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JfPtr_().oldTime();
}


Foam::surfaceScalarField& Foam::viscoBread::Jf()
{
    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::viscoBread::viscoBread
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    JPtr_(),
    JfPtr_(),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.storeOldTime();

    // Force the creation of Fs so they are read on restart
    F();
    Ff();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn
        (
            "viscoBread::viscoBread::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscoBread::~viscoBread()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscoBread::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    const volScalarField muBar(Ibar*mu_);

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*muBar*DLambda_/magSTrial));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::viscoBread::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());
    J().correctBoundaryConditions();

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // prepare DEpsilon
    const Time& time = mesh().time();
    scalar dTimeSc = time.deltaTValue();
    dimensionedScalar dTime("dTime", dimTime, dTimeSc); // -- timestep

    // -- temperature
    volScalarField T = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 273;

    // -- composition
    volScalarField alphaS = mesh().lookupObject<volScalarField>("alphaS");
    volScalarField alphaL = mesh().lookupObject<volScalarField>("alphaL");
    volScalarField alphaG = 1 - alphaL - alphaS;
    alphaG.correctBoundaryConditions();

    // -- relaxation time, Young modulus, Poisson ration, pre-elastic matrix factor
    // volScalarField tau = (9.0 * (2.0 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    volScalarField tau = (9.0 * (2.0 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    tau.correctBoundaryConditions();
    dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);

    // volScalarField K_my = K_ ;
    // volScalarField K_my = K_ * (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    // volScalarField K_my = K_ * (- Foam::atan(4e4 * alphaG - 4e3) / 0.1 + 16.71);

    // volScalarField E = 9 * mu_ * K_my  / (3 * K_my + mu_);
    // volScalarField nu = 0.5 * (3 * K_my - 2 * mu_) / (3 * K_my + mu_);

    dimensionedScalar preCoeff = 1 / ((1 + nu) * (1 - 2 * nu));
    // volScalarField preCoeff = 1 / ((1 + nu) * (1 - 2 * nu));


    volTensorField invF = inv(F());
    invF.correctBoundaryConditions();
    volTensorField S = J() * invF & sigma & invF.T();
    S.correctBoundaryConditions();
    volSymmTensorField D0 = 0 * sigma;
    D0.replace(symmTensor::XX, S.component(tensor::XX) - nu * S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    D0.replace(symmTensor::XY, (2 * nu + 2) * S.component(tensor::XY));
    D0.replace(symmTensor::XZ, (2 * nu + 2) * S.component(tensor::XZ));
    D0.replace(symmTensor::YY, - nu * S.component(tensor::XX) + S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    D0.replace(symmTensor::YZ, (2 * nu + 2) * S.component(tensor::YZ));
    D0.replace(symmTensor::ZZ, - nu * S.component(tensor::XX) - nu * S.component(tensor::YY) + S.component(tensor::ZZ));
    D0.correctBoundaryConditions();

    volTensorField dEpsPInit = 1 / J() / E / tau * dTime * (F() & D0 & F().T());
    dEpsPInit.correctBoundaryConditions();
    // DEpsilonP_ = 1 / E / tau * dTime * D0;
    // DEpsilonP_ = 1 / E / tau * dTime * sigma;
    DEpsilonP_.replace(symmTensor::XX, dEpsPInit.component(tensor::XX));
    DEpsilonP_.replace(symmTensor::XY, dEpsPInit.component(tensor::XY));
    DEpsilonP_.replace(symmTensor::XZ, dEpsPInit.component(tensor::XZ));
    DEpsilonP_.replace(symmTensor::YY, dEpsPInit.component(tensor::YY));
    DEpsilonP_.replace(symmTensor::YZ, dEpsPInit.component(tensor::YZ));
    DEpsilonP_.replace(symmTensor::ZZ, dEpsPInit.component(tensor::ZZ));
    DEpsilonP_.correctBoundaryConditions();

    const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");
    // const volVectorField& DD = mesh().lookupObject<volVectorField>("DD");
    // volTensorField gradDD = fvc::grad(DD);
    
    volSymmTensorField dEpsilon = symm(gradDD);
    dEpsilon.correctBoundaryConditions();
    volTensorField dEpsInit = 1 / J() * F() & symm(gradDD) & F().T();
    dEpsilon.replace(symmTensor::XX, dEpsInit.component(tensor::XX));
    dEpsilon.replace(symmTensor::XY, dEpsInit.component(tensor::XY));
    dEpsilon.replace(symmTensor::XZ, dEpsInit.component(tensor::XZ));
    dEpsilon.replace(symmTensor::YY, dEpsInit.component(tensor::YY));
    dEpsilon.replace(symmTensor::YZ, dEpsInit.component(tensor::YZ));
    dEpsilon.replace(symmTensor::ZZ, dEpsInit.component(tensor::ZZ));

    dEpsilon.correctBoundaryConditions();
    volSymmTensorField dSigma = E * dEpsilon;

    dSigma.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    dSigma.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XY));
    dSigma.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XZ));
    dSigma.replace(symmTensor::YY, E * preCoeff * ((nu) * dEpsilon.component(symmTensor::XX) + (1 - nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    dSigma.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::YZ));
    dSigma.replace(symmTensor::ZZ, E * preCoeff * ((nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (1 - nu) * dEpsilon.component(symmTensor::ZZ)));
    dSigma.correctBoundaryConditions();

    volSymmTensorField dSigmaP = E * DEpsilonP_;

    dSigmaP.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * DEpsilonP_.component(symmTensor::XX) + (nu) * DEpsilonP_.component(symmTensor::YY) + (nu) * DEpsilonP_.component(symmTensor::ZZ)));
    dSigmaP.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::XY));
    dSigmaP.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::XZ));
    dSigmaP.replace(symmTensor::YY, E * preCoeff * ((nu) * DEpsilonP_.component(symmTensor::XX) + (1 - nu) * DEpsilonP_.component(symmTensor::YY) + (nu) * DEpsilonP_.component(symmTensor::ZZ)));
    dSigmaP.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::YZ));
    dSigmaP.replace(symmTensor::ZZ, E * preCoeff * ((nu) * DEpsilonP_.component(symmTensor::XX) + (nu) * DEpsilonP_.component(symmTensor::YY) + (1 - nu) * DEpsilonP_.component(symmTensor::ZZ)));
    dSigmaP.correctBoundaryConditions();

    sigma = sigma.oldTime() + (dSigma - dSigmaP);
    sigma.correctBoundaryConditions();
}


void Foam::viscoBread::correct(surfaceSymmTensorField& sigma)
{
    Info << "Not implemented for surface tensor." <<endl;
}


Foam::scalar Foam::viscoBread::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().time().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        return
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonPf_.primitiveField()
                  - DEpsilonPf_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonPf_.internalField()
                  - DEpsilonPf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
#endif
    }
    else
    {
        return
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonP_.primitiveField()
                  - DEpsilonP_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
#endif
    }
}


void Foam::viscoBread::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;

    epsilonP_ += DEpsilonP_;
    epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    reduce(numCellsYielding, sumOp<int>());

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::viscoBread::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient: already done by updateF
    // if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    // {
    //     F() = fvc::average(relFf()) & F().oldTime();
    // }
    // else
    // {
    //     F() = relF() & F().oldTime();
    // }

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon(0.5*log(symm(F().T() & F())));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

    // Take reference to internal fields
#ifdef OPENFOAM_NOT_EXTEND
    const symmTensorField& DEpsilonPI = DEpsilonP_.primitiveField();
    const symmTensorField& plasticNI = plasticN_.primitiveField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().primitiveField();
    const scalarField& epsilonEqI = epsilonEq.primitiveField();
#else
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();
#endif

    // Calculate error field
    const symmTensorField DEpsilonPErrorI
    (
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL)
    );

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::viscoBread::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is lover 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


void Foam::viscoBread::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    J().writeOpt() = IOobject::AUTO_WRITE;
    bEbar_.writeOpt() = IOobject::AUTO_WRITE;

    Ff().writeOpt() = IOobject::AUTO_WRITE;
    Jf().writeOpt() = IOobject::AUTO_WRITE;
    bEbarf_.writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
