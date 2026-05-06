// -- gas conservation
fvScalarMatrix pGEq
(
    // -- accumulation
    alphaG / univR / T * fvm::ddt(Mg, J, pG)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG) + alphaG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)

    
    == 

    // -- Darcy diffusion
    - fvm::laplacian(J * Finv & jGTilda, pG)

    // -- CO2 generation 
    + J * mCO2

    // -- water evaporation
    + J * evRate

    - fvm::Sp(J * kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * yV * molMV, pG)
    + J * kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * yV * molMV * pG
);

pGEq.relax();

SolverPerformance<scalar>::debug = 1;

pResidual = pGEq.solve().initialResidual();

SolverPerformance<scalar>::debug = 0;

pG.relax();

pG.correctBoundaryConditions();

// ************************************************************************* //