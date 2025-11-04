// -- debug level -- 3 debug levels, see ifs below
int debug = 1;
pimple.dict().readIfPresent("debug", debug);

// -- writing of fields according to debug level
Info << "Debug level: " << debug <<endl;

IOobject::writeOption writeOption1 = IOobject::NO_WRITE;
IOobject::writeOption writeOption2 = IOobject::NO_WRITE;
if (debug > 0)
{
    writeOption1 = IOobject::AUTO_WRITE;
}
if (debug == 2)
{
    writeOption2 = IOobject::AUTO_WRITE;
}

#include "createFieldsFromMesh.H"

#include "createEffProp.H"

// -- gas molar mass
volScalarField Mg 
(
    IOobject
    (
        "Mg",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1 / (omegaV / molMV + (1 - omegaV) / molMC)
);

// -- water vapor molar fraction
volScalarField yV
(
    IOobject
    (
        "yV",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption1
    ),
    omegaV / molMV * Mg
);

// -- generation CO2 rate
volScalarField mCO2
(
    IOobject
    (
        "mCO2",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    R0 * Foam::exp(-(T - Tm) / deltaT)
);

// -- deformarion gradient
volTensorField F 
(
    IOobject
    (
        "F",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption1
    ),
    physics->gradD().T() + I
);

volTensorField Finv = inv(F);

// -- Jacobian
volScalarField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption1
    ),
    Foam::det(F)
);

// -- density of the gas in pores
volScalarField rhoG
(
    IOobject
    (
        "rhoG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption1
    ),
    pG * Mg / univR / T
);

// -- water vapor pressure
volScalarField pV
(
    IOobject
    (
        "pV",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption2
    ),
    pG * yV
);

// -- saturated pressure
volScalarField pSat
(
    IOobject
    (
        "pSat",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13))
);

// -- moisture content kgLwater/kGS
volScalarField moisture
(
    IOobject
    (
        "moisture",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaL * rhoL / (alphaS * rhoS)
);

// -- equalibrium water vapor pressure
volScalarField pVEq
(
    IOobject
    (
        "pVEq",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    pSat
);

// -- saturated pressure
volScalarField TSat
(
    IOobject
    (
        "TSat",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    dummyT * (-2.4432439345097229911e7+2.30650e5*Foam::log(pV / dummyP))/(-1.1598177639491068526e5+5000*Foam::log( pV / dummyP))
);

// -- evaporation rate
volScalarField evRate
(
    IOobject
    (
        "evRate",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    kMPCOpen * (1 - alphaS) * molMC / univR / T * (pVEq - pV) * 0
);

// -- permeability
volScalarField kGForPerm
(
    IOobject
    (
        "kGForPerm",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    (1 - 1.1 * (alphaL / (1 - alphaS))),
    "zeroGradient"
);

volScalarField nuG
(
    IOobject
    (
        "nuG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1.52451e-6 * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113),
    "zeroGradient"
);

// volScalarField kGForPerm = (1 - 1.1 * (alphaL / (1 - alphaS)));
// volScalarField nuG = 1.52451e-6 * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113);

// -- permeability
volScalarField permGLViscG
(
    IOobject
    (
        "permGLViscG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    permGLViscGSc / nuG * kGForPerm,
    "zeroGradient"
);

// -- Darcy gas flux
volTensorField jGTilda
(
    IOobject
    (
        "jGTilda",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    - rhoG * permGLViscG * Finv.T()
);

// -- water vapor diffusive flux
volTensorField jDVTilda
(
    IOobject
    (
        "jDVTilda",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    - rhoG * DEff * Finv.T()
);

surfaceScalarField gradPg
(
    IOobject
    (
        "gradPg",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    fvc::snGrad(pG)
);


// -- debug field of the vapor flow
volVectorField vaporVel
(
    IOobject
    (
        "vaporVel",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    - omegaV * permGLViscG * (Finv.T() & fvc::grad(pG)) - DEff * (Finv.T() & fvc::grad(yV))
);


// -- basic info about initialization of the fields 
Info << "Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
Info << "Min (pSat): " << min(pSat).value() << ", max (pSat): " << max(pSat).value() << "." << endl;
Info << "Min (rhoG): " << min(rhoG).value() << ", max (rhoG): " << max(rhoG).value() << "." << endl;
Info << "Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
Info << "univR: " << univR.value() << endl;
Info << "Min (moisture): " << min(moisture).value() << ", max (moisture): " << max(moisture) << "." << endl;
Info << "Min (permGLViscG): " << min(permGLViscG).value() << ", max (permGLViscG): " << max(permGLViscG) << "." << endl;

// ************************************************************************* //