// volVectorField preKoeff1 = - J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T())));
// ((omegaV * jDVTilda & fvc::grad(omegaV)) + (omegaC * jDCTilda & fvc::grad(omegaC)) + (omegaAir * jDATilda & fvc::grad(omegaAir)));

// -- conservation of the water wapor
fvScalarMatrix omegaVEq
(
    // -- accumulation
    // alphaG * fvm::ddt(J, rhoG, omegaV)
    // + fvm::Sp(J * rhoG * fvc::ddt(alphaG), omegaV)

    alphaG / univR / T * pG * fvm::ddt(Mg, J, omegaV)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, pG) + alphaG * pG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), omegaV)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaV, pG)
    // - fvm::div(fvc::interpolate(J * (Finv & jGTilda & fvc::grad(pG))) & physics->mesh().Sf(), omegaV)
    // - fvc::div(fvc::interpolate(J * (Finv & jGTilda & fvc::grad(pG))) & physics->mesh().Sf(), omegaV)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaV)
    // - fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg, Mg)

    // + fvc::div(J * Finv & jC)
    // + fvc::div(fvc::interpolate(J * Finv & jC) & physics->mesh().Sf(), omegaV)
    // + fvc::laplacian(J * omegaV * (Finv & jDVTilda), pV)

    // -- stable alternative to jC
    // + fvc::laplacian(J * omegaV * (Finv & jDVTilda), omegaV)
    // + fvc::laplacian(J * omegaC * (Finv & jDCTilda), omegaC)
    // + fvc::laplacian(J * omegaAir * (Finv & jDATilda), omegaAir)

    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDVTilda.prevIter()), omegaV.prevIter())
    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDCTilda.prevIter()), omegaC.prevIter())
    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDATilda.prevIter()), omegaAir.prevIter())

    // -- this I comment for 
    // + fvm::laplacian(J * omegaV * (Finv & jDVTilda), omegaV)
    // + fvc::laplacian(J * omegaV * (Finv & jDCTilda), omegaC)
    // + fvc::laplacian(J * omegaV * (Finv & jDATilda), omegaAir)

    // - fvc::div(fvc::interpolate(J * omegaV) * jPgCorr)


    - fvc::div(fvc::interpolate(J) * jLCorr)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, pV)
    // + fvc::div(fvc::interpolate(omegaV * J) * jLCorr)

    // + fvc::div(fvc::interpolate(omegaV * J) * jLCorr)

    // + fvc::div(fvc::interpolate(omegaV * J) * jLCorr)

        // - fvc::div(J * (Finv & jL))


    // + fvc::laplacian(J * Finv & jDVTilda * (Foam::pow(omegaV, 2) + Foam::pow(omegaC, 2) + Foam::pow(omegaAir, 2)) / Mg, Mg)

    // + fvm::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * pG * Mg / molMV, omegaV)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * omegaV / molMV, pG * Mg)

    // + fvm::laplacian(J * Finv * koeff / molMV * pG * Mg, omegaV)
    // + fvc::laplacian(J * Finv * koeff / molMV * omegaV, pG * Mg)

    // + fvc::laplacian(J * Finv * koeff, pV)    
    // + fvc::laplacian(J * Finv * (1 - alpha) * alphaG * Dvc * (1 - delta) / univR / T * molMV * pG * Mg / molMV, omegaV)    
    // + fvc::laplacian(J * Finv * (1 - alpha) * alphaG * Dvc * (1 - delta) / univR / T * molMV * omegaV / molMV, pG * Mg)    


    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc * aw / univR / T * molMV, pSat)

    //  - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV))))
    // - fvm::laplacian((1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() * pG * Mg / molMV ))), omegaV)
    // - fvc::laplacian((1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() * omegaV / molMV))), pG * Mg)

    // - fvm::laplacian(- J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T()))) * pG * Mg / molMV, omegaV)
    // - fvc::laplacian(- J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * dPSatdT )) * Finv.T()))) * yV, pG)

    // - fvm::laplacian(J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHT / lambdaS * dPSatdT )) * Finv.T()))) * pG, omegaV)
    // - fvc::laplacian(J * (Finv & (( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHT / lambdaS * dPSatdT )) * Finv.T()))) * omegaV, pG)
    
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