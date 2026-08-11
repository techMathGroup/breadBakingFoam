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

// -- water vapor mass fraction
// volScalarField omegaC
// (
//     IOobject
//     (
//         "omegaC",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     physics->mesh()
// );

volScalarField omegaC
(
    IOobject
    (
        "omegaC",
        runTime.timeName(),
        physics->mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0 - omegaV,
    "calculated"
);

// -- solid mass fraction
volScalarField omegaS
(
    IOobject
    (
        "omegaS",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// -- dough matrix density
volScalarField rhoD
(
    IOobject
    (
        "rhoD",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// -- gas pressure
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

// -- temperature
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

// -- crust of bread
volScalarField alphaCrust
(
    IOobject
    (
        "alphaCrust",
        runTime.timeName(),
        physics->mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    physics->mesh()
);

// ************************************************************************* //