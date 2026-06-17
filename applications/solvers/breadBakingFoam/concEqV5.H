// volVectorField preKoeff1 = - J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T())));

// -- conservation of the water wapor
fvScalarMatrix omegaVEq
(
    // -- accumulation
    alphaG * fvm::ddt(J, rhoG, omegaV)
    + fvm::Sp(J * rhoG * fvc::ddt(alphaG), omegaV)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(alpha * J * Finv & jGTilda * omegaV, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaV)
    - fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg, Mg)

    //  - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV))))
    // - fvm::laplacian((1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() * pG * Mg / molMV ))), omegaV)
    // - fvc::laplacian((1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() * omegaV / molMV))), pG * Mg)

    // - fvm::laplacian(- J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T()))) * pG * Mg / molMV, omegaV)
    // - fvc::laplacian(- J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T()))) * yV, pG)

    // - fvm::laplacian(J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHT / lambdaS * dPSatdT )) * Finv.T()))) * pG, omegaV)
    // - fvc::laplacian(J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHT / lambdaS * dPSatdT )) * Finv.T()))) * omegaV, pG)
    
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