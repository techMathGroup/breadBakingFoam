// // -- conservation of the liquid water
// fvScalarMatrix concEqL
// (
//     // -- accumulation
//     rhoL * fvm::ddt(J, alphaL)

//     ==

//     // -- evaporation source
//     - J * evRate
// );

// concEqL.relax();

// phiLResidual = concEqL.solve().initialResidual();

// alphaL.correctBoundaryConditions();

// alphaL.relax();

// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    fvm::ddt(moisture)
    // + J * alphaD * omegaS * fvc::ddt(rhoD)

    ==

    // -- evaporation source
    - J * evRate / rhoS / alphaD0
);

concEqL.relax();

phiLResidual = concEqL.solve().initialResidual();

moisture.correctBoundaryConditions();

moisture.relax();

alphaL = (moisture * alphaD * rhoS) / (rhoL + rhoS);
alphaL.correctBoundaryConditions();
alphaS = alphaD - alphaL;

// alphaL = (1-omegaS) * alphaD;
// alphaS = omegaS * alphaD;

// alphaS = (omegaS / rhoS / (omegaS / rhoS + (1 - omegaS) / rhoL)) * alphaD;
// alphaL = alphaD - alphaS;



// ************************************************************************* //