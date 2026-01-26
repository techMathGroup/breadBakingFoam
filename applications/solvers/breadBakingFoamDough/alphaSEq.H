// -- conservation of the solid in bread
// fvScalarMatrix alphaSEq
// (
//     // -- accumulation
//     rhoS * fvm::ddt(J, alphaS)
// );

// alphaSEq.relax();

// alphaSEq.solve();

// alphaS.correctBoundaryConditions();

// alphaS.relax();

fvScalarMatrix alphaSEq
(
    // -- accumulation
    // rhoS * fvm::ddt(J, alphaS)
    fvm::ddt(J, rhoD, alphaD)

    == 

    - J * evRate
);

alphaSEq.relax();

alphaSEq.solve();

alphaD.correctBoundaryConditions();

alphaD.relax();

// ************************************************************************* //