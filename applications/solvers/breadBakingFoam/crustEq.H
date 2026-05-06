volScalarField kCrust = 52*dimensionedScalar("kCrust", dimless/dimTime, 42e3) * Foam::exp(- 64151 / (8.314 * T/dummyT));


// -- crust formation
fvScalarMatrix crustEq
(
    // -- accumulation
    fvm::ddt(J,alphaCrust)

    ==

    // -- evaporation source
    J*kCrust * (1-alphaCrust)
);

crustEq.relax();

crustEq.solve().initialResidual();

alphaCrust.correctBoundaryConditions();

alphaCrust.relax();

// ************************************************************************* //