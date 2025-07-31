{
    fvScalarMatrix EEqn
    (
        // -- accumulation
        // J * (alphaS * rhoS * CpSVSF + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + (1 - omegaV) * CpCO2)) * fvm::ddt(T)

        // CpS * rhoS * fvm::ddt(J, alphaS, T)
        // + CpL * rhoL * fvm::ddt(J, alphaL, T)
        // + thermo.Cp() * rhoG * fvm::ddt(J, alphaG, T)
        // + thermo.Cp() * J * alphaG * T * fvc::ddt(rhoG)

        ==

        // -- difusion
        // - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(yV)) * CpVapor +  (1 - omegaV) * thermo.Cp() * (jGTilda & fvc::grad(pG)) - molMC / Mg * rhoG * DEff * thermo.Cp() * fvc::grad((1 - yV)) )) + 

        // -- mass convection 
        // - J * thermo.Cp() * (Finv & jGTilda & fvc::grad(pG) & fvc::grad(T)) 

        // // -- mass difusion 
        // + (J * DEff * rhoG / Mg * (Finv & fvc::grad(yV) & fvc::grad(T)) * (molMV * CpVapor - molMC * thermo.Cp()))

        // -- convection 
        // fvm::div(fvc::interpolate(thermo.Cp() * J * Finv & jGTilda & fvc::grad(pG)) & physics->mesh().Sf(), T) + 

        // -- heat conduction
        fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

        // -- evaporation
        // + J * evRate * dHConst
        + J * evRate * dHT

        // + J * mCO2 * dHS
    );

    EEqn.relax();

    TResidual = EEqn.solve().initialResidual();
}

// ************************************************************************* //