// // 2. Pomocná proměnná Theta pro plyn
// volScalarField ThetaG = alphaG * Mg / T / univR;

// volVectorField Um = physics->U();
// surfaceScalarField phiM = fvc::flux(Um);   // Tok rychlosti matrice
// surfaceScalarField phiThetaG = fvc::flux(ThetaG * Um);
// volScalarField sourceCoeff
// (
//     fvc::div(phiThetaG)
// );

// -- gas conservation
fvScalarMatrix pGEq
(
    // -- accumulation
    // alphaG / univR / T * omegaAir * fvm::ddt(Mg, J, pG)
    // + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG, omegaAir) + alphaG * omegaAir * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)

    alphaG / univR / T * fvm::ddt(Mg, J, pG)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG) + alphaG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)

    // + fvm::div(phiThetaG, pG) - fvm::Sp(fvc::div(phiThetaG), pG)
    // + fvm::Sp(sourceCoeff, pG)
    
    == 

    // -- Darcy diffusion
    // - fvm::laplacian(J * omegaAir * Finv & jGTilda, pG)
    - fvm::laplacian(J * Finv & jGTilda, pG)
    // - fvc::laplacian(J * Finv & jDCTilda / Mg, Mg)

    // - fvc::laplacian(J * Finv & jDVTilda, omegaV)
    // - fvc::laplacian(J * Finv & jDCTilda, omegaC)
    // - fvc::laplacian(J * Finv & jDATilda, omegaAir)

    // + fvc::laplacian(3 * J * omegaV * (Finv & jDVTilda), omegaV)
    // + fvc::laplacian(3 * J * omegaC * (Finv & jDCTilda), omegaC)
    // + fvc::laplacian(3 * J * omegaAir * (Finv & jDATilda), omegaAir)

    // + fvc::laplacian(J.prevIter() * (Finv.prevIter() & jDVTilda.prevIter()), omegaV.prevIter())
    // + fvc::laplacian(J.prevIter() * (Finv.prevIter() & jDCTilda.prevIter()), omegaC.prevIter())
    // + fvc::laplacian(J.prevIter() * (Finv.prevIter() & jDATilda.prevIter()), omegaAir.prevIter())

    // - fvc::laplacian(J * omegaV * (Finv & jDVTilda), omegaV)
    // - fvc::laplacian(J * omegaV * (Finv & jDCTilda), omegaC)
    // - fvc::laplacian(J * omegaV * (Finv & jDATilda), omegaAir)

    // + fvc::laplacian(J * (Finv & jDVTilda), omegaV)
    // + fvc::laplacian(J * (Finv & jDCTilda), omegaC)
    // + fvc::laplacian(J * (Finv & jDATilda), omegaAir)

    // - fvc::div(fvc::interpolate(J) * jPgCorr)

    // - fvc::div(fvc::interpolate(J) * jLCorr)
    // + fvm::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * omegaV * Mg / molMV, pG)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * pG / molMV, omegaV * Mg)
    // - fvc::div(fvc::interpolate(J) * jLCorr)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, pV)

    // + fvc::div(fvc::interpolate(J) * jLCorr)

    - fvc::div(fvc::interpolate(J) * jLCorr)

    // + fvc::div(fvc::interpolate(3 * omegaV * J) * jLCorr)

    // + fvc::div(fvc::interpolate(J * Finv & jC) & physics->mesh().Sf())

    // -- CO2 generation 
    + J * mCO2
    // + mCO2

    // -- water evaporation
    + J * evRate

    // - fvc::div(fvc::interpolate(J) * jLCorr)

    // + fvm::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * omegaV * Mg / molMV, pG)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * pG / molMV, omegaV * Mg)

    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc * aw / univR / T * molMV, pSat)

        // + fvm::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc / univR / T * molMV * yV, pG)
        // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc / univR / T * molMV * pG, yV)
        // + fvm::laplacian(J * Finv * koeff * yV, pG)
        // + fvc::laplacian(J * Finv * koeff * pG, yV)
    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc * (1 - delta) / univR / T * molMV, pV)
    // + fvc::laplacian(J * Finv * (1 - alpha) * alphaG * Dvc * (1 - delta) / univR / T * molMV, pV)     


    - fvm::Sp(J * kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * yV * molMV, pG)
    + J * kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * yV * molMV * pG
);

// pGEq.relax(pRelax);
pGEq.relax();

SolverPerformance<scalar>::debug = 1;

pResidual = pGEq.solve().initialResidual();

SolverPerformance<scalar>::debug = 0;

pG.relax();
// pG.relax(pRelax);

pG.correctBoundaryConditions();

// ************************************************************************* //