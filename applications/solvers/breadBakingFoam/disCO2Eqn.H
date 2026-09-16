// -- conservation of the dissolved CO2
fvScalarMatrix CO2inLEqn
(
    // -- accumulation
    alphaD0 * rhoSScalar * fvm::ddt(moisture, XCO2InL)

    ==

    // -- dissolution
    - fvm::Sp(J * kCO2VSF * (1e5 / kH / molMC.value() * dummyP), XCO2InL)
    + J * kCO2VSF * yC * pG
);

CO2inLEqn.relax();

CO2inLEqn.solve().initialResidual();

XCO2InL.correctBoundaryConditions();

XCO2InL.relax();

// ************************************************************************* //