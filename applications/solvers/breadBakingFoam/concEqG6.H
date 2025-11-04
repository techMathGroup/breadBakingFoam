// -- gas conservation
fvScalarMatrix pGEq
(
    // -- accumulation
    (1 - alphaS - alphaL) / univR / T * fvm::ddt(Mg, J, pG)
    // + J * pG * Mg / univR * ((1 / T) * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)) + (1 - alphaS - alphaL) * (- Foam::pow(T, - 2)) * fvc::ddt(T))
    + fvm::Sp(J * Mg / univR * ((1 / T) * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)) + (1 - alphaS - alphaL) * (- Foam::pow(T, - 2)) * fvc::ddt(T)), pG)
    
    == 

    // -- Darcy diffusion
    - fvm::laplacian(J * Finv & jGTilda, pG)

    // -- CO2 generation 
    + J * mCO2

    // -- water evaporation
    + J * evRate
);

pGEq.relax();

pResidual = pGEq.solve().initialResidual();

pG.relax();

pG.correctBoundaryConditions();

// -- Darcy gas mass flux intensity 
jGTilda = - rhoG * permGLViscG * Finv.T();

// ************************************************************************* //