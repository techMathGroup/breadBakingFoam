{
    // volScalarField h1 = (J * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir));
    // volScalarField h2 = (J * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * rhoG * alphaG;
    // volScalarField h3 = rhoG * J * alphaG * omegaV;
    // volScalarField h4 = omegaL * rhoD;




    fvScalarMatrix EEqn
    (
        (alphaD0 * rhoD * (omegaS * CpS + omegaL * CpL) + J * alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T) 
        + (dHConst + (CpVapor - CpL) * (T - TRef)) * (J * evRate)
        // + (dHConst + CpL * (T - TRef)) * (fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, pV))
        // + (dHConst + CpL * (T - TRef)) * fvc::div(J * alphaG * (1 - alpha) * Dvc * molMV / univR / T * (Finv & (Finv.T() & fvc::grad(pV))))
        // + (dHConst + CpL * (T - TRef)) * (fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, awPSat))
        // - (dHConst + CpL * (T - TRef)) * fvc::div(fvc::interpolate(J) * jLCorr)
        // - (dHConst + CpL * (T - TRef)) * fvc::div(fvc::interpolate(J) * jLCorr)
        // - (CpVapor * (T - TRef)) * (fvc::div(J * Finv & jVE))
        // - (CpVapor * (T - TRef)) * (fvc::div(fvc::interpolate(J) * jLCorr))
        // - (CpCO2 * (T - TRef)) * (fvc::div(J * Finv & jCE))
        // - (CpAir * (T - TRef)) * (fvc::div(J * Finv & jAE))
        
        // - fvc::laplacian(J * (dHConst + CpL * (T - TRef)) * (Finv & Finv.T()) * dKoeffLucas, pV)
        // + fvc::div(fvc::interpolate((dHConst + CpL * (T - TRef)) * J) * jLCorr)
        // + fvc::div(fvc::interpolate(( dHConst + CpL * (T - TRef)) * J) * jLCorr)
        // - fvc::laplacian(J * (dHConst + CpL * (T - TRef)) * (Finv & Finv.T()) * alphaG * (1 - alpha) * Dvc * molMV / univR / T, pV)
        // - fvc::laplacian(J * (dHConst + CpL * (T - TRef)) * (Finv & Finv.T()) * dKoeffLucas, awPSat)
        // + fvc::div((CpVapor * (T - TRef)) * (J * Finv & jVE))
        // + fvc::div(fvc::interpolate((CpVapor * (T - TRef)) * (J)) * jLCorr)
        // + fvc::div((CpCO2 * (T - TRef)) * (J * Finv & jCE))
        // + fvc::div((CpAir * (T - TRef)) * (J * Finv & jAE))

        == 

        fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
        // + fvc::laplacian(dKoeffLucas * (dHConst) * J * Finv & Finv.T(), awPSat)

        // alphaD0 * CpS * fvm::ddt(omegaS, rhoD, T) + 
        // alphaD0 * CpL * fvm::ddt(omegaL, rhoD, T) +
        // h1 * fvm::ddt(rhoG, alphaG, T) + 
        // (rhoG * alphaG * T) * fvc::ddt(h1) 

        // - alphaD0 * CpS * TRef * fvc::ddt(omegaS, rhoD)
        // - alphaD0 * CpL * TRef * fvc::ddt(omegaL, rhoD)
        // - TRef * fvc::ddt(h2)

        // // + dHConst * (fvc::ddt(h3) - alphaD0 * fvc::ddt(h4))
        // + dHConst * (fvc::ddt(h3))

        // // - fvm::div(fvc::interpolate(CpL * dKoeffLucas * (J * Finv & Finv.T() & fvc::grad(pV))) & physics->mesh().Sf(), T)
        // // + fvc::div(CpL * TRef * dKoeffLucas * (J * Finv & Finv.T() & fvc::grad(pV)))
        // // - fvc::div(dHConst * dKoeffLucas * (J * Finv & Finv.T() & fvc::grad(pV)))
        // // - fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * (dHConst + CpL * (T - TRef)), pV)
        // // - fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas * dHConst, pV)


        // + fvc::div(dHConst * J * Finv & jVE)
        // + fvm::div(fvc::interpolate(CpVapor * J * Finv & jVE) & physics->mesh().Sf(), T)
        // - fvc::div(CpVapor * J * Finv & jVE * TRef)

        // + fvm::div(fvc::interpolate(CpCO2 * J * Finv & jCE) & physics->mesh().Sf(), T)
        // - fvc::div(CpCO2 * J * Finv & jCE * TRef)

        // + fvm::div(fvc::interpolate(CpAir * J * Finv & jAE) & physics->mesh().Sf(), T)
        // - fvc::div(CpAir * J * Finv & jAE * TRef)


        // // -- bilance ktere muzu verit 
        // // J * (alphaS * rhoS * CpS + alphaL * rhoL * CpL + (1 - alphaS - alphaL) * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        // // J * (alphaD * rhoD * (omegaS * CpS + omegaL * CpL) + alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        // // (alphaD0 * rhoD * (omegaS * CpS + omegaL * CpL) + J * alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T)
        // // alphaD0 * rhoD * CpRef * fvm::ddt(T)


        // ==

        // // - (fvc::grad(T) & (omegaV * CpVapor * (jGTilda & fvc::grad(pG))))
        // // - (fvc::grad(T) & (CpVapor *  (jDVTilda & fvc::grad(omegaV))))   
        // // + (fvc::grad(T) & (CpVapor *  (omegaV * jC)))   
        // // // - (fvc::grad(T) & (CpVapor * omegaV / Mg * (jDVTilda & fvc::grad(Mg))))   
        
        // // - (fvc::grad(T) & (omegaC * CpCO2 * (jGTilda & fvc::grad(pG))))
        // // - (fvc::grad(T) & (CpCO2 *  (jDCTilda & fvc::grad(omegaC))))
        // // + (fvc::grad(T) & (CpCO2 *  (omegaC * jC)))   
        // // // - (fvc::grad(T) & (CpCO2 * omegaC / Mg * (jDVTilda & fvc::grad(Mg))))   

        // // - (fvc::grad(T) & (omegaAir * CpAir * (jGTilda & fvc::grad(pG))))
        // // - (fvc::grad(T) & (CpAir *  (jDATilda & fvc::grad(omegaAir))))
        // // + (fvc::grad(T) & (CpAir *  (omegaAir * jC)))   
        // // // - (fvc::grad(T) & (CpAir * omegaAir / Mg * (jDVTilda & fvc::grad(Mg))))  

        // // + (fvc::grad(T) & (CpL * dKoeffLucas * (Finv.T() & fvc::grad(pV))))
        
        // // // - (fvc::grad(T) & (CpVapor *  (1 - alpha) * J * (Finv & (((alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV )) * Finv.T() ))) & fvc::grad(pV)))  

        // // // -- heat conduction
        // // + fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
        // fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
        // // + J * evRate * dHT
        // // // + J * mCO2 * dHS

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
    // sumJiHi += fvc::interpolate((jVE * (CpVapor * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi += (jLCorr * fvc::interpolate(((CpVapor * (T - TRef))))) / mag(physics->mesh().Sf());
    // // sumJiHi += (jLCorrIntoLiquid * fvc::interpolate((-dHConst + (CpVapor * (T - TRef))))) / mag(physics->mesh().Sf());
    // sumJiHi += fvc::interpolate((jCE * (CpCO2 * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi += fvc::interpolate((jAE * (CpAir * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi -= fvc::interpolate((dKoeffLucas * (Finv.T() & fvc::grad(pV)) * (dHConst + CpL * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi -= jLCorr * fvc::interpolate((dHConst + CpL * (T - TRef))) / mag(physics->mesh().Sf());
    // sumJiHi -= jLCorr * fvc::interpolate((CpL * (T - TRef))) / mag(physics->mesh().Sf());
    // sumJiHi -= fvc::interpolate((alphaG * (1 - alpha) * Dvc * molMV / univR / T * (Finv.T() & fvc::grad(pV)) * (dHConst + CpL * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi -= fvc::interpolate((dKoeffLucas * (Finv.T() & fvc::grad(awPSat)) * (dHConst + CpL * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());
    // sumJiHi = 0 * sumJiHi;
    // sumJiHi -= fvc::interpolate((dKoeffLucas * (Finv.T() & fvc::grad(pV)) * (dHConst + CpL * (T - TRef)))) & physics->mesh().Sf() / mag(physics->mesh().Sf());

    // sumJiHi.correctBoundaryConditions();

    EEqn.relax();

    TResidual = EEqn.solve().initialResidual();

    T.relax();
    T.correctBoundaryConditions();
}

// ************************************************************************* //

