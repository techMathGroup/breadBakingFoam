// -- placeholder for effective heat conductivity
volScalarField lambdaEff
(
    IOobject
    (
        "lambdaEff",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        writeOption2
    ),
    lambdaS * omegaV,
    "calculated"
);
lambdaEff.correctBoundaryConditions();

// -- placeholder for efective diffusion coefficient
volScalarField DEffvM
(
    IOobject
    (
        "DEffvM", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    DG * omegaV,
    "calculated"
);
DEffvM.correctBoundaryConditions();

volScalarField DEffcM 
(
    IOobject
    (
        "DEffcM", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    DG * omegaV,
    "calculated"
);
DEffcM.correctBoundaryConditions();

volScalarField DEffaM 
(
    IOobject
    (
        "DEffaM", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    DG * omegaV,
    "calculated"
);
DEffaM.correctBoundaryConditions();

volScalarField Dvc 
(
    IOobject
    (
        "Dvc", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    dimensionedScalar("Dvc", DG.dimensions(), 2.15e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75),
    "calculated"
);
Dvc.correctBoundaryConditions();

volScalarField Dva 
(
    IOobject
    (
        "Dva", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("Dva", DG.dimensions(), 2.65e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75),
    "calculated"
);
Dva.correctBoundaryConditions();

volScalarField Dca 
(
    IOobject
    (
        "Dca", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("Dca", DG.dimensions(), 1.67e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75),
    "calculated"
);
Dca.correctBoundaryConditions();

volScalarField a
(
    IOobject
    (
        "a",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    physics->mesh(),
    // Foam::exp(evCoef1 * T / dummyT + evCoef2),
    dimensionedScalar("a", evCoef2.dimensions(), 1),
    "calculated"
);
a.correctBoundaryConditions();

volScalarField aw
(
    IOobject
    (
        "aw",
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    physics->mesh(),
    // Foam::exp(evCoef1 * T / dummyT + evCoef2),
    dimensionedScalar("aw", dimless, 0),
    "calculated"
);
aw.correctBoundaryConditions();



// ************************************************************************* //