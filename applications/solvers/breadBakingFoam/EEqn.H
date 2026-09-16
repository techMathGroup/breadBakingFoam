// -- enthalpy conservation
{
    fvScalarMatrix EEqn
    (
        (alphaD0 * rhoSScalar * (CpS + moisture * CpL) + J * alphaG * rhoG * (omegaV * CpVapor + omegaC * CpCO2 + omegaAir * CpAir)) * fvm::ddt(T) 
        + (dHConst + (CpVapor - CpL) * (T - TRef)) * (J * evRate)
        + (dHConst + CpL * (T - TRef)) * liquidOutFlux

        == 

        fvm::laplacian(lambdaEff * J * Finv & Finv.T(), T)
    );



    EEqn.relax();

    TResidual = EEqn.solve().initialResidual();

    T.relax();
    T.correctBoundaryConditions();
}

// -- flux of liquid and species through boundary (neglected)
// sumJiHi = 0 * sumJiHi;
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

// ************************************************************************* //

