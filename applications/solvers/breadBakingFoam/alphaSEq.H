// -- conservation of the solid in bread
fvScalarMatrix alphaSEq
(
    // -- accumulation
    rhoS * fvm::ddt(J, alphaS)
);

alphaSEq.relax();

alphaSEq.solve();

alphaS.correctBoundaryConditions();

alphaS.relax();

// ************************************************************************* //