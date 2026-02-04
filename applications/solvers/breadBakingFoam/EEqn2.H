{
    fvScalarMatrix EEqn
    (
        // -- bilance ktere muzu verit 
        J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + (1 - omegaV) * CpCO2)) * fvm::ddt(T)

        // -- accumulation
        // J * (alphaS * rhoS * CpSVSF + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
        // Product-ddt accumulation terms (implicit in T)
        // CpS * rhoS * fvm::ddt(J, alphaS, T)
        // + CpL * rhoL * fvm::ddt(J, alphaL, T)
        // + (omegaV * CpVapor * (1 - alphaS - alphaL) * fvm::ddt(J, rhoG, T))
        // + J * rhoG * T * CpVapor * (fvc::ddt(omegaV) * (1 - alphaS - alphaL) + omegaV * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)))

        // + ((1 - omegaV) * CpCO2 * (1 - alphaS - alphaL) * fvm::ddt(J, rhoG, T)) 
        // + J * rhoG * T * CpCO2 * ( - fvc::ddt(omegaV) * (1 - alphaS - alphaL) + (1 - omegaV) * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)))


        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + (1 - omegaV) * CpCO2)) * fvm::ddt(T)

        // CpS * rhoS * fvm::ddt(J, alphaS, T)
        // + CpL * rhoL * fvm::ddt(J, alphaL, T)
        // - CpCO2 * rhoG * fvm::ddt(J, alphaL, T)
        // - CpCO2 * rhoG * fvm::ddt(J, alphaS, T)
        // + CpCO2 * J * (1 - alphaL - alphaS) * T * fvc::ddt(rhoG)

        ==

        - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpVapor *  (jDVTilda & fvc::grad(omegaV))))   
        
        - (fvc::grad(T) & ((1 - omegaV) * CpCO2 * (jGTilda & fvc::grad(pG))))
        + (fvc::grad(T) & (CpCO2 *  (jDVTilda & fvc::grad(omegaV))))
        // -- difusion
        // - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(yV)) * CpVapor +  (1 - omegaV) * thermo.Cp() * (jGTilda & fvc::grad(pG)) - molMC / Mg * rhoG * DEff * thermo.Cp() * fvc::grad((1 - yV)) )) + 

        // -- mass convection 
        // - J * thermo.Cp() * (Finv & jGTilda & fvc::grad(pG) & fvc::grad(T)) 

        // // -- mass difusion 
        // + (J * DEff * rhoG / Mg * (Finv & fvc::grad(yV) & fvc::grad(T)) * (molMV * CpVapor - molMC * thermo.Cp()))

        // -- convection 
        // fvm::div(fvc::interpolate(thermo.Cp() * J * Finv & jGTilda & fvc::grad(pG)) & physics->mesh().Sf(), T) + 

        // -- Darcy convective enthalpy transport (gas mixture)
        // - fvc::div(J * (Finv & qConvDarcy))

        // // -- net diffusive enthalpy transport from counter-diffusion (H2O vs CO2)
        // - fvc::div(J * (Finv & qDiffNet))

        // -- heat conduction
        + fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

        // -- evaporation
        // - J * evRate * (dHConst + (T - dummyT * 373) * (CpVapor - CpL))
        // - J * evRate * (dHConst + (dummyT * 273) * (CpL - CpVapor))
        + J * evRate * dHT

        // + J * mCO2 * ( (T) * (CpCO2 - CpS))
    );

    EEqn.relax();

    TResidual = EEqn.solve().initialResidual();
}

// ************************************************************************* //