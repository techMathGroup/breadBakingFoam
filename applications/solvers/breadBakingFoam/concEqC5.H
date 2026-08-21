// -- conservation of the CO2
fvScalarMatrix omegaCEq
(
    // -- accumulation
    // alphaG * fvm::ddt(J, rhoG, omegaC)
    // + fvm::Sp(J * rhoG * fvc::ddt(alphaG), omegaC)

    alphaG / univR / T * pG * fvm::ddt(Mg, J, omegaC)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, pG) + alphaG * pG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), omegaC)
    // + (J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, pG) + alphaG * pG * (- Foam::pow(T, - 2)) * fvc::ddt(T)) * omegaC)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaC, pG)
    // - fvc::laplacian(J * Finv & jGTilda * omegaC, pG)
    // - fvm::div(fvc::interpolate(J * (Finv & jGTilda & fvc::grad(pG))) & physics->mesh().Sf(), omegaC)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDCTilda, omegaC)
    // - fvc::laplacian(J * Finv & jDCTilda * omegaC / Mg, Mg)

    // + fvc::div(J * Finv & jC)
    // - fvc::laplacian(J * Finv & jDVTilda * omegaC / Mg, Mg)
    // + fvm::div(fvc::interpolate(J * Finv & jC) & physics->mesh().Sf(), omegaC)
    // + fvc::laplacian(J * omegaV * (Finv & jDVTilda), omegaV)
    // + fvc::laplacian(J * omegaC * (Finv & jDCTilda), omegaC)
    // + fvc::laplacian(J * omegaAir * (Finv & jDATilda), omegaAir)

    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDVTilda.prevIter()), omegaV.prevIter())
    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDCTilda.prevIter()), omegaC.prevIter())
    // + fvc::laplacian(J.prevIter() * omegaV.prevIter() * (Finv.prevIter() & jDATilda.prevIter()), omegaAir.prevIter())

    + fvc::laplacian(J * omegaC * (Finv & jDVTilda), omegaV)
    + fvm::laplacian(J * omegaC * (Finv & jDCTilda), omegaC)
    + fvc::laplacian(J * omegaC * (Finv & jDATilda), omegaAir)

    // - fvc::div(fvc::interpolate(J * omegaC) * jPgCorr)

    // + fvc::div(fvc::interpolate(omegaC * J) * jLCorr)

    // + fvc::laplacian(J * Finv & jDVTilda * (Foam::pow(omegaV, 2) + Foam::pow(omegaC, 2) + Foam::pow(omegaAir, 2)) / Mg, Mg)
    
    // -- CO2 generation 
    + J * mCO2
    // + mCO2

    // + mCO2Dis

    + J * (kCO2VSF * (1e5 / kH / molMC.value() * dummyP) * XCO2InL)
    - fvm::Sp(J * kCO2VSF * pG * Mg / molMC, omegaC)
    // - (kCO2 * pG * Mg / molMC * omegaC)
);

omegaCEq.relax();

omegaCResidual = mag(omegaCEq.solve().initialResidual());

omegaC.relax();

omegaC.correctBoundaryConditions();

// ************************************************************************* //