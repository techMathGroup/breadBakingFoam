// -- conservation of the solid in bread
fvScalarMatrix alphaSEq
(
    // -- accumulation
    rhoS * fvm::ddt(J, alphaS)

    ==

    - J * mCO2
);

alphaSEq.relax();

alphaSEq.solve();

alphaS.correctBoundaryConditions();

alphaS.relax();

// ************************************************************************* //