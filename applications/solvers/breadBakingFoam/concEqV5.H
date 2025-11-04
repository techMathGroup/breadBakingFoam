// -- conservation of the water wapor
fvScalarMatrix omegaVEq
(
    // -- accumulation
    (1 - alphaS - alphaL) * fvm::ddt(J, rhoG, omegaV)
    + fvm::Sp(J * rhoG * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)), omegaV)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaV, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaV)
    // - fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg, Mg)
    
    // -- evaporation
    + J * evRate
);

omegaVEq.relax();

omegaVResidual = omegaVEq.solve().initialResidual();

omegaV.relax();

omegaV.correctBoundaryConditions();

// -- difusive flux of the water vapor
jDVTilda = - rhoG * DEff * Finv.T();

// -- gas molar mass
Mg = 1 / (omegaV / molMV + (1 - omegaV) / molMC);
Mg.correctBoundaryConditions();

// -- vapor molar fraction
yV = omegaV / molMV * Mg;
yV.correctBoundaryConditions();

vaporVel = - omegaV * permGLViscG * (Finv.T() & fvc::grad(pG)) - DEff * (Finv.T() & fvc::grad(yV));

// ************************************************************************* //