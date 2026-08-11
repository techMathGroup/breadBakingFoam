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
    1 - omegaV - omegaC,
    "calculated"
);
omegaAir.correctBoundaryConditions();

// forAll(omegaAir.boundaryField(), patchI)
// {
//     scalarField& omegaAirP = omegaAir.boundaryFieldRef()[patchI];
//     scalarField& omegaCP = omegaC.boundaryFieldRef()[patchI];
//     scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];

//     forAll(omegaAirP, faceI)
//     {
//         omegaAirP[faceI] = 1 - omegaVP[faceI] - omegaCP[faceI];
//     }
// }

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
    1.0 / (omegaV / molMV + omegaC / molMC + omegaAir / molMAir),
    "calculated"
);
Mg.correctBoundaryConditions();

// Info << "molMV: " << molMV << ", molMC: " << molMC << ", molMAir: " << molMAir << endl;


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
    R0*Foam::exp(-Foam::pow((T - Tm) / deltaT, 2)) * alphaD0 * rhoD * omegaS / dummyRho,
    "calculated"
);
mCO2.correctBoundaryConditions();

volTensorField F = physics->mesh().lookupObject<volTensorField>("F");

volTensorField Finv = physics->mesh().lookupObject<volTensorField>("Finv");

// -- Jacobian
volScalarField J = physics->mesh().lookupObject<volScalarField>("J");

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
    pG * Mg / univR / T,
    "calculated"
);
rhoG.correctBoundaryConditions();

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
    omegaV / molMV * Mg,
    "calculated"
);
yV.correctBoundaryConditions();

volScalarField yC
(
    IOobject
    (
        "yC",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption2
    ),
    omegaC / molMC * Mg,
    "calculated"
);
yC.correctBoundaryConditions();

volScalarField yA
(
    IOobject
    (
        "yA",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption2
    ),
    omegaAir / molMAir * Mg,
    "calculated"
);
yA.correctBoundaryConditions();

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
    dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13)),
    "calculated"
);
pSat.correctBoundaryConditions();

// -- water vapor mass fraction
volScalarField omegaL
(
    IOobject
    (
        "omegaL",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1 - omegaS,
    "calculated"
);
omegaL.correctBoundaryConditions();

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
    omegaL / omegaS,
    "calculated"
);
moisture.correctBoundaryConditions();

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
    kMPCOpen * molMC / univR / T * (pG) * 0,
    "calculated"
);
evRate.correctBoundaryConditions();

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
        IOobject::NO_WRITE
    ),
    1.52451e-6 * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113),
    "calculated"
);
nuG.correctBoundaryConditions();

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
    "calculated"
);
permGLViscG.correctBoundaryConditions();

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
    Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.1415 + 0.5,
    "calculated"
);
alpha.correctBoundaryConditions();

volScalarField tortField 
(
    IOobject
    (
        "tortField", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    tortOpen * alpha + (1 - alpha) * tortClosed,
    "calculated"
);
tortField.correctBoundaryConditions();

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
    - rhoG * permGLViscG * Finv.T(),
    "calculated"
);
jGTilda.correctBoundaryConditions();

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
    - rhoG * DEffvM * Finv.T(),
    "calculated"
);
jDVTilda.correctBoundaryConditions();

volTensorField jDCTilda
(
    IOobject
    (
        "jDCTilda",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    - rhoG * DEffcM * Finv.T(),
    "calculated"
);
jDCTilda.correctBoundaryConditions();


volTensorField jDATilda
(
    IOobject
    (
        "jDATilda",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    - rhoG * DEffaM * Finv.T(),
    "calculated"
);
jDATilda.correctBoundaryConditions();

volVectorField jC
(
    IOobject
    (
        "jC",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    // physics->mesh(),
    (omegaV * (jDVTilda & fvc::grad(omegaV)) + (omegaC * jDCTilda & fvc::grad(omegaC)) + (omegaAir * jDATilda & fvc::grad(omegaAir))),
    "calculated"
);
jC.correctBoundaryConditions();


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

volScalarField alphaD
(
    IOobject
    (
        "alphaD",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    alphaD0 / J,
    "calculated"
);
alphaD.correctBoundaryConditions();

volScalarField alphaG
(
    IOobject
    (
        "alphaG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1 - alphaD,
    "calculated"
);
alphaG.correctBoundaryConditions();

// volScalarField preCoef = - alphaG * rhoG * permGLViscG;

// surfaceScalarField CO2Flux
// (
//     IOobject
//     (
//         "CO2Flux",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption1
//     ),
//     ((fvc::interpolate(preCoef * J) / mag(physics->mesh().Sf()))) * (fvc::interpolate(fvc::grad(pG)) & (fvc::interpolate(Finv.T()) & physics->mesh().Sf()))
// );


// volScalarField rhoS
// (
//     IOobject
//     (
//         "rhoS",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     physics->mesh(),
//     rhoSScalar,
// );


// volScalarField artP
// (
//     IOobject
//     (
//         "artP",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     1e3 * rhoS / rhoSScalar * dummyP
// );

// volScalarField yield
// (
//     IOobject
//     (
//         "yield",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     physics->mesh(),
//     dimensionedScalar("yield", dimless, 0)
// );

volScalarField awPSat
(
    IOobject
    (
        "awPSat",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    aw * pSat,
    "calculated"
);
awPSat.correctBoundaryConditions();

volScalarField availSurf
(
    IOobject
    (
        "availSurf",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    physics->mesh(),
    dimensionedScalar("availSurf", dimless, 0),
    "calculated"
);
availSurf.correctBoundaryConditions();

// volVectorField jDiffNow
// (
//     IOobject
//     (
//         "jDiffNow",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV))))
// );



volScalarField pV
(
    IOobject
    (
        "pV",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    pG * yV,
    "calculated"
);
pV.correctBoundaryConditions();

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    1 - Foam::pow((alphaG), 1.0/3.0),
    "calculated"
);
delta.correctBoundaryConditions();

volVectorField jVE
(
    IOobject
    (
        "jVE",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    omegaV * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaV)) + jC,
    "calculated"
);
jVE.correctBoundaryConditions();

volVectorField jCE
(
    IOobject
    (
        "jCE",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    omegaC * (jGTilda & fvc::grad(pG)) + (jDCTilda & fvc::grad(omegaC)) + jC,
    "calculated"
);
jCE.correctBoundaryConditions();

volVectorField jAE
(
    IOobject
    (
        "jAE",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    omegaAir * (jGTilda & fvc::grad(pG)) + (jDATilda & fvc::grad(omegaAir)) + jC,
    "calculated"
);
jAE.correctBoundaryConditions();

surfaceScalarField sumJiHi
(
    IOobject
    (
        "sumJiHi",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    0 * (fvc::interpolate((jVE * (CpVapor * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf())),
    "calculated"
);
sumJiHi.correctBoundaryConditions();

volScalarField dPWdX
(
    IOobject
    (
        "dPWdX",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    ((Foam::pow(100,(2.0/nCoef)))*(Foam::pow((moisture/a),(2.0/nCoef))) - Foam::pow(10000,(1/nCoef))*Foam::pow((moisture/a),(2.0/nCoef)) + Foam::pow(100,(1.0/nCoef))*Foam::pow((moisture/a),(1.0/nCoef))/(nCoef*moisture*Foam::pow((Foam::pow(100,(1.0/nCoef))*Foam::pow((moisture/a),(1.0/nCoef)) + 1),2))) * pSat,
    "calculated"
);
dPWdX.correctBoundaryConditions();

volScalarField dPSatdT
(
    IOobject
    (
        "dPSatdT",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    3.81644e6*Foam::exp((81.443*T/dummyT - 22839.16559)/(5.*T/dummyT - 230.65))/(Foam::pow(T/dummyT - 46.13, 2)) * dummyP / dummyT,
    "calculated"
);
dPSatdT.correctBoundaryConditions();

volScalarField dKoeffLucas
(
    IOobject
    (
        "dKoeffLucas",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    (1 - alpha) * alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / Dvc) + delta * (1.0 / (rhoL * Dl) * dPWdX + dHConst / lambdaS * dPSatdT) ),
    "calculated"
);
dKoeffLucas.correctBoundaryConditions();

volScalarField kMPC
(
    IOobject
    (
        "kMPC",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha * kMPCOpen + (1 - alpha) * kMPCClosed,
    "calculated"
);
kMPC.correctBoundaryConditions();


volScalarField Sl
(
    IOobject
    (
        "Sl",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    alphaG*0,
    "calculated"
);
Sl.correctBoundaryConditions();

volVectorField jL
(
    IOobject
    (
        "jL",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    - dKoeffLucas * ((Finv.T() & fvc::grad(pV)))
);
jL.correctBoundaryConditions();

surfaceScalarField jLCorr
(
    IOobject
    (
        "jLCorr",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(Finv & jL) & physics->mesh().Sf()
);

surfaceScalarField jLCorrIntoLiquid = 0*jLCorr;



// -- basic info about initialization of the fields 
Info << "Min (pG): " << min(pG).value() << ", max (pG): " << max(pG).value() << "." << endl;
Info << "Min (pSat): " << min(pSat).value() << ", max (pSat): " << max(pSat).value() << "." << endl;
Info << "Min (rhoG): " << min(rhoG).value() << ", max (rhoG): " << max(rhoG).value() << "." << endl;
Info << "Min (T): " << min(T).value() << ", max (T): " << max(T).value() << "." << endl;
Info << "univR: " << univR.value() << endl;
Info << "Min (moisture): " << min(moisture).value() << ", max (moisture): " << max(moisture) << "." << endl;
Info << "Min (permGLViscG): " << min(permGLViscG).value() << ", max (permGLViscG): " << max(permGLViscG) << "." << endl;

// ************************************************************************* //