// -- conservation of the CO2
fvScalarMatrix omegaCEq
(
    // -- accumulation
    alphaG * fvm::ddt(J, rhoG, omegaC)
    + fvm::Sp(J * rhoG * fvc::ddt(alphaG), omegaC)

    == 

    // -- Darcy diffusion
    - fvc::laplacian(J * Finv & jGTilda * omegaC, pG)

    // -- diffusion
    - fvm::laplacian(J * Finv & jDVTilda, omegaC)
    + fvc::laplacian(J * Finv & jDVTilda * omegaC / Mg, Mg)
    
    // -- CO2 generation 
    + J * mCO2
);

omegaCEq.relax();

omegaCResidual = mag(omegaCEq.solve().initialResidual());

omegaC.relax();

omegaC.correctBoundaryConditions();

// ************************************************************************* //