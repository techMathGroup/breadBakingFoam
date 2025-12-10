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
    D0_
    (
        IOobject
        (
            "D0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    D0f_
    (
        IOobject
        (
            "D0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    T_
    (
        IOobject
        (
            "TSolid",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    alphaG_
    (
        IOobject
        (
            "alphaG",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    tau_
    (
        IOobject
        (
            "tau",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTime, 0)
    ),
    invF_
    (
        IOobject
        (
            "invF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    invFf_
    (
        IOobject
        (
            "invF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    S_
    (
        IOobject
        (
            "S",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimPressure, tensor::zero)
    ),
    Sf_
    (
        IOobject
        (
            "S",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimPressure, tensor::zero)
    ),
    dEpsPInit_
    (
        IOobject
        (
            "dEpsPInit",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    dEpsPInitf_
    (
        IOobject
        (
            "dEpsPInit",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    dEpsInit_
    (
        IOobject
        (
            "dEpsInit",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    dEpsInitf_
    (
        IOobject
        (
            "dEpsInit",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    dSigma_
    (
        IOobject
        (
            "dSigma",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    dSigmaf_
    (
        IOobject
        (
            "dSigma",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    dSigmaP_
    (
        IOobject
        (
            "dSigmaP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    dSigmaPf_
    (
        IOobject
        (
            "dSigmaP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    )
{
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
    // const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    // const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    // const volScalarField muBar(Ibar*mu_);

    // // Magnitude of the deviatoric trial stress
    // const volScalarField magSTrial
    // (
    //     max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    // );

    // // Calculate scaling factor
    // const volScalarField scaleFactor(1.0 - (2.0*muBar/magSTrial));

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
            mesh(),
            (4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            zeroGradientFvPatchScalarField::typeName
            // scaleFactor*(4.0/3.0)*mu_ + K_
            // (4.0/3.0)*mu_ + K_
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
    // DEpsilonP_.storePrevIter();

    // prepare DEpsilon
    const Time& time = mesh().time();
    scalar dTimeSc = time.deltaTValue();
    dimensionedScalar dTime("dTime", dimTime, dTimeSc); // -- timestep

    // -- temperature
    const volScalarField& TItself = mesh().lookupObject<volScalarField>("T");
    T_ = TItself / dimensionedScalar("dummyT", dimTemperature, 1) - 273;

    // -- composition
    const volScalarField& alphaS = mesh().lookupObject<volScalarField>("alphaS");
    const volScalarField& alphaL = mesh().lookupObject<volScalarField>("alphaL");
    alphaG_ = 1 - alphaL - alphaS;

    // -- relaxation time, Young modulus, Poisson ration, pre-elastic matrix factor
    tau_ = (9.0 * (2.0 / 3.14 * Foam::atan((T_ - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(4e4 * alphaG_ - 4e3) / 1e-3 + 1571.75);
    dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);

    // volScalarField K_my = K_ ;
    // volScalarField K_my = K_ * (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    // volScalarField K_my = K_ * (- Foam::atan(4e4 * alphaG - 4e3) / 0.1 + 16.71);

    // volScalarField E = 9 * mu_ * K_my  / (3 * K_my + mu_);
    // volScalarField nu = 0.5 * (3 * K_my - 2 * mu_) / (3 * K_my + mu_);

    dimensionedScalar preCoeff = 1 / ((1 + nu) * (1 - 2 * nu));
    // volScalarField preCoeff = 1 / ((1 + nu) * (1 - 2 * nu));


    invF_ = inv(F());
    S_ = J() * invF_ & sigma & invF_.T();
    D0_.replace(symmTensor::XX, S_.component(tensor::XX) - nu * S_.component(tensor::YY) - nu * S_.component(tensor::ZZ));
    D0_.replace(symmTensor::XY, (2 * nu + 2) * S_.component(tensor::XY));
    D0_.replace(symmTensor::XZ, (2 * nu + 2) * S_.component(tensor::XZ));
    D0_.replace(symmTensor::YY, - nu * S_.component(tensor::XX) + S_.component(tensor::YY) - nu * S_.component(tensor::ZZ));
    D0_.replace(symmTensor::YZ, (2 * nu + 2) * S_.component(tensor::YZ));
    D0_.replace(symmTensor::ZZ, - nu * S_.component(tensor::XX) - nu * S_.component(tensor::YY) + S_.component(tensor::ZZ));

    dEpsPInit_ = 1 / J() / E / tau_ * dTime * (F() & D0_ & F().T());
    // DEpsilonP_ = 1 / E / tau * dTime * D0;
    // DEpsilonP_ = 1 / E / tau * dTime * sigma;
    // DEpsilonP_.replace(symmTensor::XX, dEpsPInit_.component(tensor::XX));
    // DEpsilonP_.replace(symmTensor::XY, dEpsPInit_.component(tensor::XY));
    // DEpsilonP_.replace(symmTensor::XZ, dEpsPInit_.component(tensor::XZ));
    // DEpsilonP_.replace(symmTensor::YY, dEpsPInit_.component(tensor::YY));
    // DEpsilonP_.replace(symmTensor::YZ, dEpsPInit_.component(tensor::YZ));
    // DEpsilonP_.replace(symmTensor::ZZ, dEpsPInit_.component(tensor::ZZ));

    const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");
    // const volVectorField& DD = mesh().lookupObject<volVectorField>("DD");
    // volTensorField gradDD = fvc::grad(DD);
    
    // volSymmTensorField dEpsilon = symm(gradDD);
    dEpsInit_ = 1 / J() * F() & symm(gradDD) & F().T();
    // dEpsilon_.replace(symmTensor::XX, dEpsInit_.component(tensor::XX));
    // dEpsilon_.replace(symmTensor::XY, dEpsInit_.component(tensor::XY));
    // dEpsilon_.replace(symmTensor::XZ, dEpsInit_.component(tensor::XZ));
    // dEpsilon_.replace(symmTensor::YY, dEpsInit_.component(tensor::YY));
    // dEpsilon_.replace(symmTensor::YZ, dEpsInit_.component(tensor::YZ));
    // dEpsilon_.replace(symmTensor::ZZ, dEpsInit_.component(tensor::ZZ));

    // volSymmTensorField dSigma = E * dEpsilon;

    dSigma_.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * dEpsInit_.component(tensor::XX) + (nu) * dEpsInit_.component(tensor::YY) + (nu) * dEpsInit_.component(tensor::ZZ)));
    dSigma_.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInit_.component(tensor::XY));
    dSigma_.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInit_.component(tensor::XZ));
    dSigma_.replace(symmTensor::YY, E * preCoeff * ((nu) * dEpsInit_.component(tensor::XX) + (1 - nu) * dEpsInit_.component(tensor::YY) + (nu) * dEpsInit_.component(tensor::ZZ)));
    dSigma_.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInit_.component(tensor::YZ));
    dSigma_.replace(symmTensor::ZZ, E * preCoeff * ((nu) * dEpsInit_.component(tensor::XX) + (nu) * dEpsInit_.component(tensor::YY) + (1 - nu) * dEpsInit_.component(tensor::ZZ)));

    // volSymmTensorField dSigmaP = E * DEpsilonP_;

    dSigmaP_.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * dEpsPInit_.component(tensor::XX) + (nu) * dEpsPInit_.component(tensor::YY) + (nu) * dEpsPInit_.component(tensor::ZZ)));
    dSigmaP_.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInit_.component(tensor::XY));
    dSigmaP_.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInit_.component(tensor::XZ));
    dSigmaP_.replace(symmTensor::YY, E * preCoeff * ((nu) * dEpsPInit_.component(tensor::XX) + (1 - nu) * dEpsPInit_.component(tensor::YY) + (nu) * dEpsPInit_.component(tensor::ZZ)));
    dSigmaP_.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInit_.component(tensor::YZ));
    dSigmaP_.replace(symmTensor::ZZ, E * preCoeff * ((nu) * dEpsPInit_.component(tensor::XX) + (nu) * dEpsPInit_.component(tensor::YY) + (1 - nu) * dEpsPInit_.component(tensor::ZZ)));

    sigma = sigma.oldTime() + (dSigma_ - dSigmaP_);
}


void Foam::viscoBread::correct(surfaceSymmTensorField& sigma)
{
        // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    Jf() = det(Ff());
    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    // DEpsilonP_.storePrevIter();

    // prepare DEpsilon
    const Time& time = mesh().time();
    scalar dTimeSc = time.deltaTValue();
    dimensionedScalar dTime("dTime", dimTime, dTimeSc); // -- timestep

    // -- temperature
    const volScalarField& TItself = mesh().lookupObject<volScalarField>("T");
    T_ = TItself / dimensionedScalar("dummyT", dimTemperature, 1) - 273;

    // -- composition
    const volScalarField& alphaS = mesh().lookupObject<volScalarField>("alphaS");
    const volScalarField& alphaL = mesh().lookupObject<volScalarField>("alphaL");
    alphaG_ = 1 - alphaL - alphaS;

    // -- relaxation time, Young modulus, Poisson ration, pre-elastic matrix factor
    tau_ = (9.0 * (2.0 / 3.14 * Foam::atan((T_ - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(4e4 * alphaG_ - 4e3) / 1e-3 + 1571.75);
    dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);
    dimensionedScalar preCoeff = 1 / ((1 + nu) * (1 - 2 * nu));

    invFf_ = inv(Ff());
    Sf_ = Jf() * invFf_ & sigma & invFf_.T();
    D0f_.replace(symmTensor::XX, Sf_.component(tensor::XX) - nu * Sf_.component(tensor::YY) - nu * Sf_.component(tensor::ZZ));
    D0f_.replace(symmTensor::XY, (2 * nu + 2) * Sf_.component(tensor::XY));
    D0f_.replace(symmTensor::XZ, (2 * nu + 2) * Sf_.component(tensor::XZ));
    D0f_.replace(symmTensor::YY, - nu * Sf_.component(tensor::XX) + Sf_.component(tensor::YY) - nu * Sf_.component(tensor::ZZ));
    D0f_.replace(symmTensor::YZ, (2 * nu + 2) * Sf_.component(tensor::YZ));
    D0f_.replace(symmTensor::ZZ, - nu * Sf_.component(tensor::XX) - nu * Sf_.component(tensor::YY) + Sf_.component(tensor::ZZ));

    dEpsPInitf_ = 1 / Jf() / E / fvc::interpolate(tau_) * dTime * (Ff() & D0f_ & Ff().T());

    const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");
    // const surfaceTensorField& gradD = mesh().lookupObject<surfaceTensorField>("grad(D)f");
    // surfaceTensorField gradDDf = gradD - gradD.oldTime();
    surfaceTensorField gradDDf = fvc::interpolate(gradDD);
    // const volVectorField& DD = mesh().lookupObject<volVectorField>("DD");
    // volTensorField gradDD = fvc::grad(DD);
    
    // const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");
    // const volVectorField& DD = mesh().lookupObject<volVectorField>("DD");
    // volTensorField gradDD = fvc::grad(DD);
    
    // volSymmTensorField dEpsilon = symm(gradDD);
    dEpsInitf_ = 1 / Jf() * Ff() & symm(gradDDf) & Ff().T();
    // dEpsilon_.replace(symmTensor::XX, dEpsInit_.component(tensor::XX));
    // dEpsilon_.replace(symmTensor::XY, dEpsInit_.component(tensor::XY));
    // dEpsilon_.replace(symmTensor::XZ, dEpsInit_.component(tensor::XZ));
    // dEpsilon_.replace(symmTensor::YY, dEpsInit_.component(tensor::YY));
    // dEpsilon_.replace(symmTensor::YZ, dEpsInit_.component(tensor::YZ));
    // dEpsilon_.replace(symmTensor::ZZ, dEpsInit_.component(tensor::ZZ));

    // volSymmTensorField dSigma = E * dEpsilon;

    dSigmaf_.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * dEpsInitf_.component(tensor::XX) + (nu) * dEpsInitf_.component(tensor::YY) + (nu) * dEpsInitf_.component(tensor::ZZ)));
    dSigmaf_.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInitf_.component(tensor::XY));
    dSigmaf_.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInitf_.component(tensor::XZ));
    dSigmaf_.replace(symmTensor::YY, E * preCoeff * ((nu) * dEpsInitf_.component(tensor::XX) + (1 - nu) * dEpsInitf_.component(tensor::YY) + (nu) * dEpsInitf_.component(tensor::ZZ)));
    dSigmaf_.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsInitf_.component(tensor::YZ));
    dSigmaf_.replace(symmTensor::ZZ, E * preCoeff * ((nu) * dEpsInitf_.component(tensor::XX) + (nu) * dEpsInitf_.component(tensor::YY) + (1 - nu) * dEpsInitf_.component(tensor::ZZ)));

    // volSymmTensorField dSigmaP = E * DEpsilonP_;

    dSigmaPf_.replace(symmTensor::XX, E * preCoeff * ((1 - nu) * dEpsPInitf_.component(tensor::XX) + (nu) * dEpsPInitf_.component(tensor::YY) + (nu) * dEpsPInitf_.component(tensor::ZZ)));
    dSigmaPf_.replace(symmTensor::XY, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInitf_.component(tensor::XY));
    dSigmaPf_.replace(symmTensor::XZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInitf_.component(tensor::XZ));
    dSigmaPf_.replace(symmTensor::YY, E * preCoeff * ((nu) * dEpsPInitf_.component(tensor::XX) + (1 - nu) * dEpsPInitf_.component(tensor::YY) + (nu) * dEpsPInitf_.component(tensor::ZZ)));
    dSigmaPf_.replace(symmTensor::YZ, E * preCoeff * (1 - 2 * nu) / 2 * dEpsPInitf_.component(tensor::YZ));
    dSigmaPf_.replace(symmTensor::ZZ, E * preCoeff * ((nu) * dEpsPInitf_.component(tensor::XX) + (nu) * dEpsPInitf_.component(tensor::YY) + (1 - nu) * dEpsPInitf_.component(tensor::ZZ)));

    sigma = sigma.oldTime() + (dSigmaf_ - dSigmaPf_);
}


Foam::scalar Foam::viscoBread::residual()
{
    Info << "Residual not implemented" << endl;
    return 1;
//     // Calculate residual based on change in plastic strain increment
//     if
//     (
//         mesh().time().lookupObject<fvMesh>
//         (
//             baseMeshRegionName()
//         ).foundObject<surfaceTensorField>("Ff")
//     )
//     {
//         return
// #ifdef OPENFOAM_NOT_EXTEND
//             gMax
//             (
//                 mag
//                 (
//                     DEpsilonPf_.primitiveField()
//                   - DEpsilonPf_.prevIter().primitiveField()
//                 )
//             )/gMax(SMALL + mag(DEpsilonPf_.prevIter().primitiveField()));
// #else
//             gMax
//             (
//                 mag
//                 (
//                     DEpsilonPf_.internalField()
//                   - DEpsilonPf_.prevIter().internalField()
//                 )
//             )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
// #endif
//     }
//     else
//     {
//         return
// #ifdef OPENFOAM_NOT_EXTEND
//             gMax
//             (
//                 mag
//                 (
//                     DEpsilonP_.primitiveField()
//                   - DEpsilonP_.prevIter().primitiveField()
//                 )
//             )/gMax(SMALL + mag(DEpsilonP_.prevIter().primitiveField()));
// #else
//             gMax
//             (
//                 mag
//                 (
//                     DEpsilonP_.internalField()
//                   - DEpsilonP_.prevIter().internalField()
//                 )
//             )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
// #endif
//     }
}


void Foam::viscoBread::updateTotalFields()
{
    // Info<< nl << "Updating total accumulated fields" << endl;

    // epsilonP_ += DEpsilonP_;
    // epsilonPf_ += DEpsilonPf_;

    // // Count cells actively yielding
    // int numCellsYielding = 0;

    // reduce(numCellsYielding, sumOp<int>());

    // Info<< "    " << numCellsYielding << " cells are actively yielding"
    //     << nl << endl;
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
//     const volSymmTensorField epsilon(0.5*log(symm(F().T() & F())));

//     // Calculate equivalent strain, for normalisation of the error
//     const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

//     // Take reference to internal fields
// #ifdef OPENFOAM_NOT_EXTEND
//     const symmTensorField& DEpsilonPI = DEpsilonP_.primitiveField();
//     const scalarField& epsilonEqI = epsilonEq.primitiveField();
// #else
//     const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
//     const scalarField& epsilonEqI = epsilonEq.internalField();
// #endif

//     // Calculate error field
//     const symmTensorField DEpsilonPErrorI
//     (
//         Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
//        /(epsilonEqI + SMALL)
//     );

//     // Max error
//     const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

//     if (maxMagDEpsilonPErr > SMALL)
//     {
//         Info<< "    " << name() << ": max time integration error = "
//             << maxMagDEpsilonPErr
//             << endl;

//         if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
//         {
//             WarningIn
//             (
//                 "Foam::scalar Foam::viscoBread::newDeltaT()"
//                 " const"
//             )   << "The error in the plastic strain is lover 50 times larger "
//                 << "than the desired value!\n    Consider starting the "
//                 << "simulation with a smaller initial time-step" << endl;
//         }

//         // Calculate the time-step scaling factor, where maxDeltaErr_ is the
//         // maximum allowed error
//         const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

//         // Return the new time-step size
//         return scaleFac*mesh().time().deltaTValue();
//     }

    return mesh().time().endTime().value();
}


void Foam::viscoBread::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    J().writeOpt() = IOobject::AUTO_WRITE;

    Ff().writeOpt() = IOobject::AUTO_WRITE;
    Jf().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
