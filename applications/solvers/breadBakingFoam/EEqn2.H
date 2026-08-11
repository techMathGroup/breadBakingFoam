{
    // volVectorField jCpTotal
    // (
    //     CpVapor   * (omegaV   * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaV))   + (omegaV / Mg)   * (jDVTilda & fvc::grad(Mg)))
    //     + CpCO2     * (omegaC   * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaC))   + (omegaC / Mg)   * (jDVTilda & fvc::grad(Mg)))
    //     + CpAir     * (omegaAir * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaAir)) + (omegaAir / Mg) * (jDVTilda & fvc::grad(Mg)))
    // );

    // // // 2. Interpolujeme vektor na stěny buněk (surfaceScalarField) pro konvekční schémata
    // surfaceScalarField phiCp = fvc::flux(jCpTotal);

    // // 3. Samotná energetická rovnice
    // fvScalarMatrix EEqn
    // (
    //     // -- Akumulace (v referenčních souřadnicích)
    //     J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)

    //     // -- Implicitní entalpický transport (přesunuto z RHS na LHS)
    //     // Používáme nekonzervativní vyjádření transportu (u . grad(T)) zapsané pomocí OpenFOAM identit
    //     + fvm::div(phiCp, T)
    //     - fvm::Sp(fvc::div(phiCp), T)

    //     ==

    //     // -- Vedení tepla (Tvůj skvěle transformovaný Total Lagrangian Laplacian)
    //      fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

    //     // -- Latentní teplo odparu (opraveno znaménko na mínus pro chlazení)
    //     + J * evRate * dHT
    // );

    // volScalarField J_rhoCpEff
    // (
    //     "J_rhoCpEff",
    //     J * (
    //         alphaS * rhoS * CpS 
    //         + alphaL * rhoL * CpL 
    //         + alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)
    //     )
    // );

    // fvScalarMatrix EEqn
    // (
    //     // -- PLNĚ IMPLICITNÍ konzervativní akumulace tepla
    //     fvm::ddt(J_rhoCpEff, T)

    //     // -- Implicitní entalpický transport (konvekce plynů)
    //     // + fvm::div(phiCp, T)
    //     // - fvm::Sp(fvc::div(phiCp), T)

    //     ==

    //     // -- Vedení tepla (Total Lagrangian Laplacian)
    //     fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

    //     - fvm::div(CpVapor * ((fvc::interpolate(J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV)))))) & physics->mesh().Sf()), T)
    //     - fvm::div(CpCO2 * ((fvc::interpolate(J * (Finv & (omegaC * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaC)))))) & physics->mesh().Sf()), T)
    //     - fvm::div(CpAir * ((fvc::interpolate(J * (Finv & (omegaAir * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaAir)))))) & physics->mesh().Sf()), T)

    //     // + fvc::div(fvc::interpolate(dHT * jDiffNew) & physics->mesh().Sf())

    //     // -- Latentní teplo odparu
    //     + J * evRate * dHT
    // );

    // dimensionedScalar Tref("Tref", dimTemperature, 273.15);
    // volScalarField CpRhoVapor = rhoG * alphaG * J * omegaV;
    // volScalarField CpRhoCO2 = rhoG * alphaG * J * omegaC;

    // volVectorField jDiffNowCO2 = -(((J * (Finv & (omegaC * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaC)))))) );
    // volVectorField jDiffNowAir = -(((J * (Finv & (omegaAir * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaAir)))))) );

    // Info << "min fvc::ddt(CpRhoVapor)" << min(fvc::ddt(CpRhoVapor)).value() << endl;
    // Info << "max fvc::ddt(CpRhoVapor)" << max(fvc::ddt(CpRhoVapor)).value() << endl;

    // fvScalarMatrix EEqn
    // (
    //     // -- PLNĚ IMPLICITNÍ konzervativní akumulace tepla
    //     // fvm::ddt(J_rhoCpEff, T)

    //     // -- Implicitní entalpický transport (konvekce plynů)
    //     // + fvm::div(phiCp, T)
    //     // - fvm::Sp(fvc::div(phiCp), T)


    //     J_rhoCpEff * fvm::ddt(T)

    //     + (T - Tref) * fvc::ddt(J_rhoCpEff)

    //     + dHConst * fvc::ddt(CpRhoVapor)

    //     + dHS * fvc::ddt(CpRhoCO2)

    //     ==

    //     // -- Vedení tepla (Total Lagrangian Laplacian)
    //     fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

    //     - fvm::div(CpVapor * fvc::interpolate(jDiffNow) & physics->mesh().Sf(), T)
    //     + fvc::div(CpVapor * Tref * fvc::interpolate(jDiffNew) & physics->mesh().Sf())

    //     - fvm::div(CpCO2 * fvc::interpolate(jDiffNowCO2) & physics->mesh().Sf(), T)
    //     + fvc::div(CpCO2 * Tref * fvc::interpolate(jDiffNowCO2) & physics->mesh().Sf())

    //     - fvm::div(CpAir * fvc::interpolate(jDiffNowAir) & physics->mesh().Sf(), T)
    //     + fvc::div(CpAir * Tref * fvc::interpolate(jDiffNowAir) & physics->mesh().Sf())

    //     - fvc::div(fvc::interpolate(dHConst * jDiffNow) & physics->mesh().Sf())
    //     - fvc::div(fvc::interpolate(dHS * jDiffNowCO2) & physics->mesh().Sf())

    //     // -- Latentní teplo odparu
    //     // + J * evRate * dHT
    // );


    // fvScalarMatrix EEqn
    // (
    //     // -- bilance ktere muzu verit 
    //     J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)

    //     // rhoS * CpS * fvm::ddt(J, alphaS, T)
    //     // + rhoL * CpL * fvm::ddt(J, alphaL, T)
    //     // - rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir) * fvm::ddt(J, alphaS, T)
    //     // - rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir) * fvm::ddt(J, alphaL, T)
    //     // + CpVapor * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaV)
    //     // + CpCO2 * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaC)
    //     // + CpAir * J * (1 - alphaS - alphaL) * T * fvc::ddt(rhoG, omegaAir)
        

    //     // -- accumulation
    //     // J * (alphaS * rhoS * CpSVSF + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
    //     // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
    //     // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * thermo.Cp()) * fvm::ddt(T)
    //     // Product-ddt accumulation terms (implicit in T)
    //     // CpS * rhoS * fvm::ddt(J, alphaS, T)
    //     // + CpL * rhoL * fvm::ddt(J, alphaL, T)
    //     // + (omegaV * CpVapor * (1 - alphaS - alphaL) * fvm::ddt(J, rhoG, T))
    //     // + J * rhoG * T * CpVapor * (fvc::ddt(omegaV) * (1 - alphaS - alphaL) + omegaV * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)))

    //     // + ((1 - omegaV) * CpCO2 * (1 - alphaS - alphaL) * fvm::ddt(J, rhoG, T)) 
    //     // + J * rhoG * T * CpCO2 * ( - fvc::ddt(omegaV) * (1 - alphaS - alphaL) + (1 - omegaV) * (- fvc::ddt(alphaS) - fvc::ddt(alphaL)))


    //     // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + (1 - omegaV) * CpCO2)) * fvm::ddt(T)

    //     // CpS * rhoS * fvm::ddt(J, alphaS, T)
    //     // + CpL * rhoL * fvm::ddt(J, alphaL, T)
    //     // - CpCO2 * rhoG * fvm::ddt(J, alphaL, T)
    //     // - CpCO2 * rhoG * fvm::ddt(J, alphaS, T)
    //     // + CpCO2 * J * (1 - alphaL - alphaS) * T * fvc::ddt(rhoG)

    //     ==

    //     // - fvc::laplacian(J * Finv & jGTilda * omegaV * CpVapor * T, pG)
    //     // - fvc::laplacian(J * Finv & jGTilda * omegaC * CpCO2 * T, pG)
    //     // - fvc::laplacian(J * Finv & jGTilda * omegaAir * CpAir * T, pG)

    //     // // // -- diffusion
    //     // - fvc::laplacian(J * Finv & jDVTilda * CpVapor * T, omegaV)
    //     // - fvc::laplacian(J * Finv & jDVTilda * omegaV / Mg * CpVapor * T, Mg)

    //     // - fvc::laplacian(J * Finv & jDVTilda * CpCO2 * T, omegaC)
    //     // - fvc::laplacian(J * Finv & jDVTilda * omegaC / Mg * CpCO2 * T, Mg)

    //     // - fvc::laplacian(J * Finv & jDVTilda * CpAir * T, omegaAir)
    //     // - fvc::laplacian(J * Finv & jDVTilda * omegaAir / Mg * CpAir * T, Mg)

    //     - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG))))
    //     - (fvc::grad(T) & (CpVapor *  (jDVTilda & fvc::grad(omegaV))))   
    //     - (fvc::grad(T) & (CpVapor * omegaV / Mg * (jDVTilda & fvc::grad(Mg))))   
        
    //     - (fvc::grad(T) & (omegaC * CpCO2 * (jGTilda & fvc::grad(pG))))
    //     - (fvc::grad(T) & (CpCO2 *  (jDVTilda & fvc::grad(omegaC))))
    //     - (fvc::grad(T) & (CpCO2 * omegaC / Mg * (jDVTilda & fvc::grad(Mg))))   

    //     - (fvc::grad(T) & (omegaAir * CpAir * (jGTilda & fvc::grad(pG))))
    //     - (fvc::grad(T) & (CpAir *  (jDVTilda & fvc::grad(omegaAir))))
    //     - (fvc::grad(T) & (CpAir * omegaAir / Mg * (jDVTilda & fvc::grad(Mg))))   
    //     // -- difusion
    //     // - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(yV)) * CpVapor +  (1 - omegaV) * thermo.Cp() * (jGTilda & fvc::grad(pG)) - molMC / Mg * rhoG * DEff * thermo.Cp() * fvc::grad((1 - yV)) )) + 

    //     // -- mass convection 
    //     // - J * thermo.Cp() * (Finv & jGTilda & fvc::grad(pG) & fvc::grad(T)) 

    //     // // -- mass difusion 
    //     // + (J * DEff * rhoG / Mg * (Finv & fvc::grad(yV) & fvc::grad(T)) * (molMV * CpVapor - molMC * thermo.Cp()))

    //     // -- convection 
    //     // fvm::div(fvc::interpolate(thermo.Cp() * J * Finv & jGTilda & fvc::grad(pG)) & physics->mesh().Sf(), T) + 

    //     // -- Darcy convective enthalpy transport (gas mixture)
    //     // - fvc::div(J * (Finv & qConvDarcy))

    //     // // -- net diffusive enthalpy transport from counter-diffusion (H2O vs CO2)
    //     // - fvc::div(J * (Finv & qDiffNet))

    //     // -- heat conduction
    //     + fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)

    //     // -- evaporation
    //     // - J * evRate * (dHConst + (T - dummyT * 373) * (CpVapor - CpL))
    //     // - J * evRate * (dHConst + (dummyT * 273) * (CpL - CpVapor))
    //     + J * evRate * dHT

    //     // + J * mCO2 * ( (T) * (CpCO2 - CpS))
    // );

    // volScalarField CpRef = dummyCp * (1132 + 5.8 * (T / dummyT - 273) + 3700 * moisture - 2310 * moisture * moisture);
// surfaceScalarField jCCorr = 0 * fvc::interpolate(J * dKoeffLucas * (Finv & (Finv.T() & fvc::grad(pV)))) & physics->mesh().Sf();
// surfaceScalarField SfMag = physics->mesh().magSf();
// forAll(jCCorr.boundaryField(), patchI)
// {
//     scalarField& jCCorrP = jCCorr.boundaryFieldRef()[patchI];
//     scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
//     scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
//     const scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];

//     forAll(jCCorrP, faceI)
//     {
//         jCCorrP[faceI] = -0.1 * rhoGP[faceI] * (omegaVP[faceI]  - 7e-3)*magSfP[faceI];
//     }
// }


    fvScalarMatrix EEqn
    (
        // -- bilance ktere muzu verit 
        // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        // J * (alphaD * rhoD * (omegaS * CpS + omegaL * CpL) + alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        (alphaD0 * rhoD * (omegaS * CpS + omegaL * CpL) + J * alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        // alphaD0 * rhoD * CpRef * fvm::ddt(T)


        ==

        - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpVapor *  (jDVTilda & fvc::grad(omegaV))))
        // + (fvc::grad(T) & (CpVapor *  (omegaV * jC)))   
        // - (fvc::grad(T) & (CpVapor * omegaV / Mg * (jDVTilda & fvc::grad(Mg))))   
        
        - (fvc::grad(T) & (omegaC * CpCO2 * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpCO2 *  (jDCTilda & fvc::grad(omegaC))))
        // + (fvc::grad(T) & (CpCO2 *  (omegaC * jC)))   
        // - (fvc::grad(T) & (CpCO2 * omegaC / Mg * (jDVTilda & fvc::grad(Mg))))   

        - (fvc::grad(T) & (omegaAir * CpAir * (jGTilda & fvc::grad(pG))))
        - (fvc::grad(T) & (CpAir *  (jDATilda & fvc::grad(omegaAir))))

        + (fvc::grad(T) & ((CpVapor + CpAir + CpCO2) * omegaV * (jDVTilda) & fvc::grad(omegaV)))
        + (fvc::grad(T) & ((CpVapor + CpAir + CpCO2) * omegaC * (jDCTilda) & fvc::grad(omegaV)))
        + (fvc::grad(T) & ((CpVapor + CpAir + CpCO2) * omegaAir * (jDATilda) & fvc::grad(omegaAir)))   

        // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * (dHConst), pV)

        // + (fvc::grad(T) & (CpAir *  (omegaAir * jC)))   
        // - (fvc::grad(T) & (CpAir * omegaAir / Mg * (jDVTilda & fvc::grad(Mg))))  

        + (fvc::grad(T) & (CpL * dKoeffLucas * (Finv.T() & fvc::grad(pV))))

        // + fvc::div(dHConst * jCCorr)
        // + fvc::div(dHConst * dKoeffLucas * (J * Finv & Finv.T() & fvc::grad(pV)))
        // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * dHConst, pV)
        // + fvc::div(jCCorr * dHConst)
        
        // - (fvc::grad(T) & (CpVapor *  (1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() ))) & fvc::grad(pV)))  

        // -- heat conduction
        // + fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
        + fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
        + J * evRate * dHT
        // + J * mCO2 * dHS

    );

    // volTensorField jDATilda = - rhoG * DEffaM * Finv.T();

    // volVectorField jDiffNowCO2 = -((J * (Finv & (omegaC * ((jGTilda & fvc::grad(pG)) + (jDCTilda & fvc::grad(omegaC)))))) );
    // volVectorField jDiffNowVap = -((J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaV)))))) );
    // volVectorField jDiffNowAir = -(((J * (Finv & (omegaAir * ((jGTilda & fvc::grad(pG)) + (jDATilda & fvc::grad(omegaAir)))))) ));
    // jDiffNow = - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV))));

    // dimensionedScalar Tref("Tref", dimTemperature, 300);

    // // dimensionedScalar dHV("dHV", dH.dimensions(), 2300000);


    // dimensionedScalar Tref("Tref", dimTemperature, 300);
    sumJiHi = 0 * sumJiHi;
    dimensionedScalar TRef("TRef", dimTemperature, 300);


    // volVectorField jVE = omegaV * (jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(omegaV)) + jC; 
    // volVectorField jCE = omegaC * (jGTilda & fvc::grad(pG)) + (jDCTilda & fvc::grad(omegaC)) + jC; 
    // volVectorField jAE = omegaAir * (jGTilda & fvc::grad(pG)) + (jDATilda & fvc::grad(omegaAir)) + jC; 
    // sumJiHi += fvc::interpolate((jVE * (CpVapor * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi += fvc::interpolate((jCE * (CpCO2 * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi += fvc::interpolate((jAE * (CpAir * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    

    sumJiHi.correctBoundaryConditions();

    EEqn.relax();

    TResidual = EEqn.solve().initialResidual();

    T.relax();
    T.correctBoundaryConditions();
}

// ************************************************************************* //

