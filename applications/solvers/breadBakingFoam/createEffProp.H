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
    lambdaS * alphaS,
    "zeroGradient"
);

// -- placeholder for efective diffusion coefficient
volScalarField DEff 
(
    IOobject
    (
        "DEff", 
        runTime.timeName(),
        physics->mesh(),
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    DG * (1 - alphaS),
    "zeroGradient"
);

// ************************************************************************* //