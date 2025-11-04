// -- vol fraction solid
volScalarField alphaS
(
    IOobject
    (
        "alphaS",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// -- vol fraction fluid
volScalarField alphaL
(
    IOobject
    (
        "alphaL",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// -- water vapor mass fraction
volScalarField omegaV
(
    IOobject
    (
        "omegaV",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// gas pressure
volScalarField pG
(
    IOobject
    (
        "pG",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// ************************************************************************* //