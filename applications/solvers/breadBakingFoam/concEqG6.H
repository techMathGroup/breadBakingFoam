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
    alphaG / univR / T * fvm::ddt(Mg, J, pG)
    + fvm::Sp(J * Mg / univR * ((1 / T) * fvc::ddt(alphaG) + alphaG * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)

    // + fvm::div(phiThetaG, pG) - fvm::Sp(fvc::div(phiThetaG), pG)
    // + fvm::Sp(sourceCoeff, pG)
    
    == 

    // -- Darcy diffusion
    - fvm::laplacian(J * Finv & jGTilda, pG)

    - fvc::laplacian(J * Finv & jDVTilda / Mg, Mg)

    // -- CO2 generation 
    + J * mCO2

    // -- water evaporation
    + J * evRate

    - fvm::Sp(J * kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * yV * molMV, pG)
    + J * kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * yV * molMV * pG
);

// pGEq.relax(pRelax);
pGEq.relax();

// SolverPerformance<scalar>::debug = 1;

pResidual = pGEq.solve().initialResidual();

// SolverPerformance<scalar>::debug = 0;

pG.relax();
// pG.relax(pRelax);

pG.correctBoundaryConditions();

// ************************************************************************* //