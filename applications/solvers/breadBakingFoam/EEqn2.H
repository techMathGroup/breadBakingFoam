{
    fvScalarMatrix EEqn
    (
        // -- bilance ktere muzu verit 
        J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)

        // rhoS * CpS * fvm::ddt(J, alphaS, T)
        // + rhoL * CpL * fvm::ddt(J, alphaL, T)
        // - rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir) * fvm::ddt(J, alphaS, T)
        // - rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir) * fvm::ddt(J, alphaL, T)
        // + CpVapor * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaV)
        // + CpCO2 * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaC)
        // + CpAir * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaAir)
        

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

        // - fvc::laplacian(J * Finv & jGTilda * omegaV * CpVapor * T, pG)
        // - fvc::laplacian(J * Finv & jGTilda * omegaC * CpCO2 * T, pG)
        // - fvc::laplacian(J * Finv & jGTilda * omegaAir * CpAir * T, pG)

        // // // -- diffusion
        // - fvc::laplacian(J * Finv & jDVTilda * CpVapor * T, omegaV)
        // - fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg * CpVapor * T, Mg)

        // - fvc::laplacian(J * Finv & jDVTilda * CpCO2 * T, omegaC)
        // - fvc::laplacian(J * Finv & jDVTilda * omegaC / Mg * CpCO2 * T, Mg)

        // - fvc::laplacian(J * Finv & jDVTilda * CpAir * T, omegaAir)
        // - fvc::laplacian(J * Finv & jDVTilda * omegaAir / Mg * CpAir * T, Mg)

        - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpVapor *  (jDVTilda & fvc::grad(omegaV))))   
        - (fvc::grad(T) & (CpVapor * omegaV / Mg * (jDVTilda & fvc::grad(Mg))))   
        
        - (fvc::grad(T) & (omegaC * CpCO2 * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpCO2 *  (jDVTilda & fvc::grad(omegaC))))
        - (fvc::grad(T) & (CpCO2 * omegaC / Mg * (jDVTilda & fvc::grad(Mg))))   

        - (fvc::grad(T) & (omegaAir * CpAir * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpAir *  (jDVTilda & fvc::grad(omegaAir))))
        - (fvc::grad(T) & (CpAir * omegaAir / Mg * (jDVTilda & fvc::grad(Mg))))   
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