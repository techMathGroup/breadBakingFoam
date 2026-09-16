// -- conservation of the CO2
fvScalarMatrix omegaCEq
(
    // -- accumulation
    alphaG / univR / T * pG * fvm::ddt(Mg, J, omegaC)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, pG) + alphaG * pG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), omegaC)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaC, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDCTilda, omegaC)

    // -- correction on mixture-averaged diffusion
    + fvc::laplacian(J * omegaC * (Finv & jDVTilda), omegaV)
    + fvm::laplacian(J * omegaC * (Finv & jDCTilda), omegaC)
    + fvc::laplacian(J * omegaC * (Finv & jDATilda), omegaAir)
    
    // -- CO2 generation by fermentation
    + J * mCO2

    // -- CO2 dissolution
    + J * (kCO2VSF * (1e5 / kH / molMC.value() * dummyP) * XCO2InL)
    - fvm::Sp(J * kCO2VSF * pG * Mg / molMC, omegaC)
);

omegaCEq.relax();

omegaCResidual = mag(omegaCEq.solve().initialResidual());

omegaC.relax();

omegaC.correctBoundaryConditions();


// ************************************************************************* //