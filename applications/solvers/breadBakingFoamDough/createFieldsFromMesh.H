// -- vol fraction solid


// -- vol fraction fluid
// volScalarField omegaS
// (
//     IOobject
//     (
//         "omegaS",
//         runTime.timeName(),
//         physics->mesh(),
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     physics->mesh()
// );



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