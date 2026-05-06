// -- conservation of the water wapor
fvScalarMatrix omegaVEq
(
    // -- accumulation
    alphaG * fvm::ddt(J, rhoG, omegaV)
    + fvm::Sp(J * rhoG * fvc::ddt(alphaG), omegaV)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaV, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaV)
    + fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg, Mg)
    
    // -- evaporation
    + J * evRate

    - fvm::Sp(J * kMPC * Mg / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * pG, omegaV)
    + J * kMPC * Mg / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * pG * omegaV
);

omegaVEq.relax();
omegaVResidual = omegaVEq.solve().initialResidual();

omegaV.relax();

omegaV.correctBoundaryConditions();

// ************************************************************************* //