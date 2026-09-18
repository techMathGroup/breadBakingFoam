// -- gas phase conservation
fvScalarMatrix pGEq
(
    // -- accumulation
    alphaG / univR / T * fvm::ddt(Mg, J, pG)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG) + alphaG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)
    
    == 

    // -- Darcy diffusion
    - fvm::laplacian(J * Finv & jGTilda, pG)

    // -- CO2 dissolution 
    + J * (kCO2VSF * (1e5 / kH / molMC.value() * dummyP) * XCO2InL)
    - fvm::Sp(J * kCO2VSF * yC, pG)

    + fvc::laplacian(J * (Finv & jDVTilda), omegaV)
    + fvc::laplacian(J * (Finv & jDCTilda), omegaC)
    + fvc::laplacian(J * (Finv & jDATilda), omegaAir)

    // -- CO2 generation by fermentation 
    + J * mCO2

    // -- water evaporation
    + J * evRate

    // -- linearization of evaporation
    - fvm::Sp(J * kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * yV * molMV, pG)
    + J * kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * yV * molMV * pG
);

pGEq.relax();

SolverPerformance<scalar>::debug = 1;

pResidual = pGEq.solve().initialResidual();

SolverPerformance<scalar>::debug = 0;

pG.relax();

pG.correctBoundaryConditions();

// ************************************************************************* //