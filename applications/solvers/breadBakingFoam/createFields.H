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


// -- generation CO2 rate
volScalarField nCO2
(
    IOobject
    (
        "nCO2",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    mCO2 / molMC
    // mCO2
);

// -- deformarion gradient
volTensorField F 
(
    IOobject
    (
        "solverF",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    physics->gradD().T() + I
);

volTensorField Finv = inv(F);

// -- Jacobian
volScalarField J
(
    IOobject
    (
        "solverJ",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
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

// -- density of the gas in pores
volScalarField yV
(
    IOobject
    (
        "yV",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption2
    ),
    omegaV
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
    kMPCOpen * (1 - alphaS) * molMC / univR / T * (pG) * 0
);

// -- permeability
// -- legacy from Zhang
// volScalarField kGForPerm
// (
//     IOobject
//     (
//         "kGForPerm",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     (1 - 1.1 * (alphaL / (1 - alphaS))),
//     "zeroGradient"
// );

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
    permGLViscGSc / nuG,
    "zeroGradient"
);

// -- pores opening function
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.14 + 0.5
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
        IOobject::NO_WRITE
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

volScalarField omegaAir
(
    IOobject
    (
        "omegaAir",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1 - omegaV - omegaC
);


// volVectorField fluxG
// (
//     IOobject
//     (
//         "fluxG",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     - rhoG * permGLViscG * fvc::grad(pG)
// );

volScalarField alphaG
(
    IOobject
    (
        "alphaG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1 - alphaS - alphaL
);

volScalarField preCoef = - alphaG * rhoG * permGLViscG;

surfaceScalarField CO2Flux
(
    IOobject
    (
        "CO2Flux",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption1
    ),
    ((fvc::interpolate(preCoef * J) / mag(physics->mesh().Sf()))) * (fvc::interpolate(fvc::grad(pG)) & (fvc::interpolate(Finv.T()) & physics->mesh().Sf()))
);


// surfaceScalarField gradPg
// (
//     IOobject
//     (
//         "gradPg",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     fvc::snGrad(pG)
// );
// Info << "Here6" <<endl;

// -- debug field of the vapor flow
// volVectorField vaporVel
// (
//     IOobject
//     (
//         "vaporVel",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     - omegaV * permGLViscG * (Finv.T() & fvc::grad(pG)) - DEff * (Finv.T() & fvc::grad(yV))
// );


// -- basic info about initialization of the fields 
Info << "Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
Info << "Min (pSat): " << min(pSat).value() << ", max (pSat): " << max(pSat).value() << "." << endl;
Info << "Min (rhoG): " << min(rhoG).value() << ", max (rhoG): " << max(rhoG).value() << "." << endl;
Info << "Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
Info << "univR: " << univR.value() << endl;
Info << "Min (moisture): " << min(moisture).value() << ", max (moisture): " << max(moisture) << "." << endl;
Info << "Min (permGLViscG): " << min(permGLViscG).value() << ", max (permGLViscG): " << max(permGLViscG) << "." << endl;

// ************************************************************************* //