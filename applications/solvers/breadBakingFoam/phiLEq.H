// volScalarField a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
// volScalarField alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.14 + 0.5;
// volScalarField kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;

volScalarField citatel = J * kMPC / univR * molMV * Foam::pow(100, 1.0 / nCoef) * Foam::pow((alphaL*rhoL)/(alphaS*rhoS*a), 1.0 / nCoef) * pSat;
volScalarField zavorka = Foam::pow(100, 1.0 / nCoef) * Foam::pow((alphaL*rhoL)/(alphaS*rhoS*a), 1.0 / nCoef) + 1.0;

// volScalarField citatel = J * kMPC / univR * molMV * Foam::pow(100, 1.0 / nCoef) * Foam::pow(availMois, 1.0 / nCoef) * pSat;
// volScalarField zavorka = Foam::pow(100, 1.0 / nCoef) * Foam::pow(availMois, 1.0 / nCoef) + 1.0;

// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    rhoL * fvm::ddt(J, alphaL)

    ==

    // - fvc::laplacian(J * Finv & Finv.T() & jGTilda * (1 - alpha) * alphaG * DG / univR / T * molMV, pSat)

    // -- evaporation source
    - J * evRate

    - fvm::Sp(citatel / (Foam::pow(T / dummyT, 0.5) * dummyT * Foam::pow(zavorka, 2) * nCoef * alphaL), alphaL)
    + citatel / (Foam::pow(T / dummyT, 0.5) * dummyT * Foam::pow(zavorka, 2) * nCoef)
);

// SolverPerformance<scalar> solverPerfD;


concEqL.relax();

phiLResidual = concEqL.solve().initialResidual();

alphaL.correctBoundaryConditions();

alphaL.relax();

// ************************************************************************* //