// -- conservation of the water wapor
fvScalarMatrix omegaVEq
(
    // -- accumulation
    alphaG / univR / T * pG * fvm::ddt(Mg, J, omegaV)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, pG) + alphaG * pG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), omegaV)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaV, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaV)

    // -- correction on diffusion
    + fvm::laplacian(J * omegaV * (Finv & jDVTilda), omegaV)
    + fvc::laplacian(J * omegaV * (Finv & jDCTilda), omegaC)
    + fvc::laplacian(J * omegaV * (Finv & jDATilda), omegaAir)
    
    // -- evaporation
    + J * evRate

    - fvm::Sp(J * kMPC * availSurf * Mg / univR / Foam::pow(T / dummyT, 0.5) / dummyT * pG, omegaV)
    + J * kMPC * availSurf * Mg / univR / Foam::pow(T / dummyT, 0.5) / dummyT * pG * omegaV
);

omegaVEq.relax();
omegaVResidual = omegaVEq.solve().initialResidual();

omegaV.relax();

omegaV.correctBoundaryConditions();

// ************************************************************************* //