// -- conservation of the solid in bread
fvScalarMatrix alphaSEq
(
    // -- accumulation of the solid part of the dough matrix phase
    alphaD0 * fvm::ddt(rhoD, omegaS)
);

alphaSEq.relax();

alphaSEq.solve();

omegaS.correctBoundaryConditions();

omegaS.relax();

// ************************************************************************* //