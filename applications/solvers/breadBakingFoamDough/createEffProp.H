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
    lambdaS * moisture,
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
    DG * Foam::pow((1 - 1.11 * moisture), 1.333),
    // DG * Foam::pow((1 - 1.11 * ((1-omegaS) * alphaD) / ((1 - omegaS) * alphaD + (1 - alphaD))) * (1 - alphaS), 1.333),
    "zeroGradient"
);

// ************************************************************************* //