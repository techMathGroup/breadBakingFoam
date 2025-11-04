// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    rhoL * fvm::ddt(J, alphaL)

    ==

    // -- evaporation source
    - J * evRate
);

concEqL.relax();

phiLResidual = concEqL.solve().initialResidual();

alphaL.correctBoundaryConditions();

alphaL.relax();

// ************************************************************************* //