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
            "alphaG_",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    tau1VSF_
    (
        IOobject
        (
            "tau1VSF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTime, 0),
        "zeroGradient"
    ),
    tau2VSF_
    (
        IOobject
        (
            "tau2VSF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTime, 0),
        "zeroGradient"
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
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
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
    // dEpsPInit_
    // (
    //     IOobject
    //     (
    //         "dEpsPInit",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedTensor("zero", dimless, tensor::zero)
    // ),
    dSigmaTensP_
    (
        IOobject
        (
            "dSigmaTensP_",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimPressure, symmTensor::zero)
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
    dE_
    (
        IOobject
        (
            "dE",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
    SEl1_
    (
        IOobject
        (
            "SEl1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    SEl2_
    (
        IOobject
        (
            "SEl2",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    SV1_
    (
        IOobject
        (
            "SV1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    SV2_
    (
        IOobject
        (
            "SV2",
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
    ),
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
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    hard_
    (
        IOobject
        (
            "hard",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor(1,1,1,1,1,1,1,1,1))
    ),
    dS_
    (
        IOobject
        (
            "dS",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    EField_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0),
        "zeroGradient"
    ),
    nuField_
        (
        IOobject
        (
            "nu",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0),
        "zeroGradient"
    ),
    tau1_(0),
    tau2_(0),
    tGelat_(0),
    alphaBiot_(1.0),
    mu0Raw_(dimensionedScalar("mu0Raw", dimPressure, 0)),
    kappa0Raw_(dimensionedScalar("kappa0Raw", dimPressure, 0)),
    muV1Raw_(dimensionedScalar("muV1Raw", dimPressure, 0)),
    muV2Raw_(dimensionedScalar("muV2Raw", dimPressure, 0)),
    kappaVRaw_(dimensionedScalar("kappaVRaw", dimPressure, 0)),
    mu0Baked_(dimensionedScalar("mu0Baked", dimPressure, 0)),
    kappa0Baked_(dimensionedScalar("kappa0Baked", dimPressure, 0))
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

    tau1_ = dict.lookupOrDefault<scalar>("tau1", 1.0);
    tau2_ = dict.lookupOrDefault<scalar>("tau2", 1.0);
    tGelat_ = dict.lookupOrDefault<scalar>("tGelat", 65.0);

    alphaBiot_ = dict.lookupOrDefault<scalar>("alphaBiot", 1.0);

    // Info << "mu0 Raw: " << mu0Raw_ << endl;

    // --- Raw State Parameters (Ambient / Proofing / Early Baking) ---
    mu0Raw_ = dict.lookupOrDefault<dimensionedScalar>("mu0Raw", dimensionedScalar("mu0Raw", dimPressure, 300));         // 0.3 kPa
    kappa0Raw_ = dict.lookupOrDefault<dimensionedScalar>("kappa0Raw", dimensionedScalar("kappa0Raw", dimPressure, 1000));   // 1.0 kPa

    muV1Raw_ = dict.lookupOrDefault<dimensionedScalar>("muV1Raw", dimensionedScalar("muV1Raw", dimPressure, 3000));       // 3.0 kPa
    muV2Raw_ = dict.lookupOrDefault<dimensionedScalar>("muV2Raw", dimensionedScalar("muV2Raw", dimPressure, 3000));       // 3.0 kPa
    kappaVRaw_ = dict.lookupOrDefault<dimensionedScalar>("kappaVRaw", dimensionedScalar("kappaVRaw", dimPressure, 6000)); // 6.0 kPa
    // dimensionedScalar tau_raw("tau_raw", dimTime, 20);                // 20 seconds

    // --- Baked State Parameters (Post-Gelatinization T > 75°C) ---
    mu0Baked_ = dict.lookupOrDefault<dimensionedScalar>("mu0Baked", dimensionedScalar("mu0Baked", dimPressure, 12500));    // 25.0 kPa
    kappa0Baked_ = dict.lookupOrDefault<dimensionedScalar>("kappa0Baked", dimensionedScalar("kappa0Baked", dimPressure, 40000)); // 85.0 kPa

    Info << "mu0 Raw: " << mu0Raw_ << endl;
    Info << "kappa0Raw_ Raw: " << kappa0Raw_ << endl;
    Info << "muV1Raw_ Raw: " << muV1Raw_ << endl;
    Info << "muV2Raw_ Raw: " << muV2Raw_ << endl;
    Info << "Tau1: " << tau1_ << endl;
    Info << "Tau2: " << tau2_ << endl;
    Info << "TGelat: " << tGelat_ << endl;
    SEl1_.storeOldTime();
    SEl2_.storeOldTime();
    SV1_.storeOldTime();
    SV2_.storeOldTime();
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
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    J() = det(F());
    J().correctBoundaryConditions();

    volSymmTensorField CSymm = symm(T(F()) & F());
    volSymmTensorField invC = inv(CSymm);

    dimensionedScalar dTime("dTime", dimTime, mesh().time().deltaTValue()); 
    dimensionedScalar pRef("pRef", dimPressure, 1e5);   

    alphaG_ = mesh().lookupObject<volScalarField>("alphaG"); // Normalize alphaG to be dimensionless

    // -- temperature
    T_ = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 273;
    const volScalarField pG = mesh().lookupObject<volScalarField>("pG");

    // -- Compute current gelatinization state based on temperature field T
    // volScalarField kappaGel = 0.5 * (1.0 + Foam::tanh((T_ - 50) / 10.0));
    volScalarField kappaGel = 0.5 + Foam::atan((T_ - 65) / 2.0) / 3.14159; // Smooth transition from 0 to 1 around T = 65°C

    // -- Dynamically blend the elastostatic ground state (from soft paste to rigid crumb)
    volScalarField mu0 = (1-alphaG_) * ((1.0 - kappaGel) * mu0Raw_ + kappaGel * mu0Baked_);
    volScalarField kappa0 = (1-alphaG_) * (((1.0 - kappaGel) * kappa0Raw_ + kappaGel * kappa0Baked_)* (1 + (100 - 1) * (0.5 - Foam::atan((alphaG_ - 0.05) / 0.01) / 3.141592)));
    // volScalarField kappa0 = ((1.0 - kappaGel) * kappa0Raw_ + kappaGel * kappa0Baked_);

    volScalarField muV1 = (1.0 - kappaGel) * muV1Raw_ + kappaGel * muV1Raw_;
    volScalarField muV2 = (1.0 - kappaGel) * muV2Raw_ + kappaGel * muV2Raw_;
    volScalarField lambdaV = (1.0 - kappaGel) * kappaVRaw_ + kappaGel * kappaVRaw_;

    // tau_ = tau0_ / (1.0 - kappaGel + 1e-5) * dimensionedScalar("dummyTime", dimTime, 1); // Avoid division by zero
    tau1VSF_ = tau1_ * (2.0 * (1.0 - kappaGel) + 20.0 * kappaGel) * dimensionedScalar("dummyTime", dimTime, 1); // Avoid division by zero
    // scalar logTauRaw = Foam::log(2.0);
    // scalar logTauBaked = Foam::log(200.0);
    // tau1VSF_ = tau1_ * Foam::exp(logTauRaw + kappaGel * (logTauBaked - logTauRaw)) * dimensionedScalar("dummyTime", dimTime, 1);
    tau2VSF_ = tau2_ * (2.0 * (1.0 - kappaGel) + 20.0 * kappaGel) * dimensionedScalar("dummyTime", dimTime, 1); // Avoid division by zero
    tau1VSF_.correctBoundaryConditions();
    tau2VSF_.correctBoundaryConditions();

    // -- Baseline isotropic Neo-Hookean response (Contains the purely elastic bulk compression)
    // Info << "mu0: " << min(J()) << endl;
    volSymmTensorField SInf = mu0 * (I - invC) + kappa0 * Foam::log(J()) * invC;

    // Info << "mu2: " << min(J()) << endl;

    SEl1_ = muV1 * (I - invC) + lambdaV * Foam::log(J()) * invC;
    SEl1_.correctBoundaryConditions();

    // Info << "mu3: " << min(J()) << endl;

    SEl2_ = muV2 * (I - invC) + lambdaV * Foam::log(J()) * invC;
    SEl2_.correctBoundaryConditions();

    // Info << "mu4: " << min(J())  << endl;

    // Calculate pure elastic increment against the true stored old time
    volSymmTensorField deltaS0elastic1 = SEl1_ - SEl1_.oldTime();
    volSymmTensorField deltaS0elastic2 = SEl2_ - SEl2_.oldTime();

    // Integration coefficients
    volScalarField oden1 = dTime / tau1VSF_;
    volScalarField oden2 = dTime / tau2VSF_;
    volScalarField expFactor1 = Foam::exp(-oden1);
    volScalarField integrationFactor1 = (1.0 - expFactor1) / oden1;
    volScalarField expFactor2 = Foam::exp(-oden2);
    volScalarField integrationFactor2 = (1.0 - expFactor2) / oden2;

    // -- BIOT POROMECHANICAL COUPLING & TOTAL STRESS ASSEMBLY
    // volSymmTensorField S_total = SInf - alphaBiot_ * J() * (pG - pRef) * invC;
    volSymmTensorField S_total = SInf;

    // if (mesh().time().timeOutputValue() > 10)
    // {
        // Update the internal viscous state variable
        SV1_ = expFactor1 * SV1_.oldTime() + integrationFactor1 * deltaS0elastic1;
        SV1_.correctBoundaryConditions();
        S_total += SV1_ ;

        // SV2_ = expFactor2 * SV2_.oldTime() + integrationFactor2 * deltaS0elastic2;
        // SV2_.correctBoundaryConditions();
        // S_total += SV2_ ;
    // }
    S_total.correctBoundaryConditions();

    // Push-forward Second Piola-Kirchhoff (S_total) to Cauchy Spatial Stress (sigma)
    sigma = symm( (1.0 / J()) * (F() & S_total & T(F())) );
    sigma.correctBoundaryConditions();
}


void Foam::viscoBread::correct(surfaceSymmTensorField& sigma)
{
    
}


Foam::scalar Foam::viscoBread::residual()
{
    // Info << "Residual not implemented" << endl;
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
