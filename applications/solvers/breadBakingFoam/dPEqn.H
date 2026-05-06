// -- conservation of the liquid water
// fvScalarMatrix dPEqn
// (
//     // -- accumulation
//     fvm::ddt(dP)

//     ==

//     0.63662 / Foam::pow(dP, 2) * univR * T / Mg / pG * 3.14 / 6 * fvc::(rhoG * Foam::pow(dP, 3))
// );

// SolverPerformance<scalar> solverPerfD;


// dPEqn.relax();

// dPEqn.solve().initialResidual();

// dP.correctBoundaryConditions();

// dP.relax();

// ************************************************************************* //