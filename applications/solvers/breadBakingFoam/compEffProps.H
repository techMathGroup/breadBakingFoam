// -- calculation of the effective model parameters

// -- binary coefficients
Dvc = dimensionedScalar("Dvc", DG.dimensions(), 2.15e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dva = dimensionedScalar("Dva", DG.dimensions(), 2.65e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dca = dimensionedScalar("Dca", DG.dimensions(), 1.67e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dvc.correctBoundaryConditions();
Dva.correctBoundaryConditions();
Dca.correctBoundaryConditions();

// -- gas dynamic viscosity
nuG = gasDynamicViscosityPrefactor * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113); // -- gas viscosity at correct temperature
nuG.correctBoundaryConditions();

// -- effective diffusivity
DEffvM = alphaG * alpha / tortOpen * 1.0 / ((yC / Dvc + yA / Dva) + yV / (1 - omegaV) * (omegaC / Dvc + omegaAir / Dva));
DEffcM = alphaG * alpha / tortOpen * 1.0 / ((yV / Dvc + yA / Dca) + yC / (1 - omegaC) * (omegaV / Dvc + omegaAir / Dca));
DEffaM = alphaG * alpha / tortOpen * 1.0 / ((yC / Dca + yV / Dva) + yA / (1 - omegaAir) * (omegaC / Dca + omegaV / Dva));
DEffvM.correctBoundaryConditions();
DEffcM.correctBoundaryConditions();
DEffaM.correctBoundaryConditions();

// -- logarithmic interpolation constants
scalar logKRaw = Foam::log(gasPermeabilityRaw.value());
scalar logKBaked = Foam::log(gasPermeabilityBaked.value());

// volScalarField alphaGP = 0.5 + Foam::atan((T/dummyT - 273 - 65) / 10.0) / 3.14159; // Smooth transition from 0 to 1 around T = 65°C

// -- permiability MPC Zhang 2025
gasDarcyMobility = Foam::exp(logKRaw + alpha * (logKBaked - logKRaw)) * dimensionedScalar("gasPermeability", gasPermeabilityRaw.dimensions(), 1.0) / nuG;
// gasDarcyMobility = Foam::exp(logKRaw + alphaGP * (logKBaked - logKRaw)) * dimensionedScalar("gasPermeability", gasPermeabilityRaw.dimensions(), 1.0) / nuG;
gasDarcyMobility.correctBoundaryConditions();

// -- effective heat conductivity https://www.sciencedirect.com/science/article/pii/S0960308508000333
// -- ratio of wall thickness to wall-pore thickness (delta)
// lambdaEff = dummyLamdba *  (lambdaS.value() - 0.134 * alphaG - 0.219 * Foam::pow(alphaG, 0.5));

forAll(lambdaEff, i)
{
    if (alphaG[i] > SMALL)
    {
        lambdaEff[i] = lambdaS.value() - 0.134 * alphaG[i] - 0.219 * Foam::pow(alphaG[i], 0.5);
        delta[i] = 1 - Foam::pow((alphaG[i]), 1.0/3.0);
    }
    else
    {
        lambdaEff[i] = lambdaS.value();
        delta[i] = 1 - SMALL;
    }
}

forAll(lambdaEff.boundaryField(), patchI)
{
    scalarField& lambdaEffP = lambdaEff.boundaryFieldRef()[patchI];
    scalarField& alphaGP = alphaG.boundaryFieldRef()[patchI];
    scalarField& deltaP = delta.boundaryFieldRef()[patchI];
    forAll(lambdaEffP, faceI)
    {
        if (alphaGP[faceI] > SMALL)
        {
            lambdaEffP[faceI] = (lambdaS.value() - 0.134 * alphaGP[faceI] - 0.219 * Foam::pow(alphaGP[faceI], 0.5));
            deltaP[faceI] = 1 - Foam::pow((alphaGP[faceI]), 1.0/3.0);
        }
        else
        {
            lambdaEffP[faceI] = lambdaS.value();
            deltaP[faceI] = 1 - SMALL;
        }
    }
}

// -- calculation of variables for Lucas closed-cell vapor flux
volScalarField dPVdT = rhoG * univR / Mg * yV;
dPVdT.correctBoundaryConditions();

// -- derivation of saturated water pressure with respect to moisture content (Oswin activity)
dPSatAwWdX = ((Foam::pow(100,(2.0/nCoef)))*(Foam::pow((moistureSafe/a),(2.0/nCoef))) - Foam::pow(10000,(1.0/nCoef))*Foam::pow((moistureSafe/a),(2.0/nCoef)) + Foam::pow(100,(1.0/nCoef))*Foam::pow((moistureSafe/a),(1.0/nCoef))/(nCoef*moistureSafe*Foam::pow((Foam::pow(100,(1.0/nCoef))*Foam::pow((moistureSafe/a),(1.0/nCoef)) + 1),2))) * pSat;
dPSatAwWdX.correctBoundaryConditions();

// volScalarField DlVSF = dimensionedScalar("dummyDl", Dl.dimensions(), 1) * 9.75*Foam::exp(-7689.9/(T/dummyT)-0.0858*moisture);
// dKoeffLucas = (1 - alpha) * alphaG / ((1 - delta) * (univR * T / molMV  * (1 - yV) / Dvc) + delta * (1.0 / (rhoL * DlVSF) * dPSatAwWdX + dHConst / lambdaS * dPSatAwdT) );
dKoeffLucas = (1 - alpha) * alphaG / ((1 - delta) * (univR * T / molMV  * (1 - yV) / Dvc) + delta * (1.0 / (rhoL * Dl) * dPSatAwWdX + dHConst / lambdaS * dPSatAwdT) );
dKoeffLucas.correctBoundaryConditions();

// -- henrys law parameters and dissolution kinetics
dimensionedScalar KH0("KH0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.035);
volScalarField kH = KH0 * Foam::exp(2400.0 * (1.0 / (T / dummyT) - 1.0 / 298.15));
dimensionedScalar kCO2("kCO2", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.015);
volScalarField kCO2VSF = kCO2 * availSurf / 8.314 / Foam::pow(T / dummyT, 0.5) * molMC.value();
mCO2Dis = kCO2VSF * (1e5 * XCO2InL / kH / molMC.value() * dummyP - yC * pG);
pCSat = 1e5 * XCO2InL / kH / molMC.value() * dummyP;

// -- purlis lambda alternative (legacy)
// lambdaEff.correctBoundaryConditions();
// lambdaEff = dummyLamdba *  (lambdaS.value());
// lambdaEff = dummyLamdba *  (lambdaS.value() * alphaD + 0.02 * alphaG);

// -- lucas lambda alternative (legacy)
// volScalarField phi = (moisture - 0.11111) / (0.8518 - 0.11111);
// lambdaEff = dummyLamdba *  (0.423 * phi + 0.18 * (1-phi));
// lambdaEff = dummyLamdba *  (0.37 * phi + 0.18 * (1-phi));
// lambdaEff = dummyLamdba *  0.3;

// -- derivation of saturated water pressure with respect to temperature (lostie activity -- legacy)
// dPSatAwdT = 3.81644e6*Foam::exp((81.443*T/dummyT - 22839.16559)/(5.*T/dummyT - 230.65))/(Foam::pow(T/dummyT - 46.13, 2)) * dummyP / dummyT * aw;
// dPSatAwdT.correctBoundaryConditions();
// dPSatAwWdX = 0.06230065326*Foam::exp(-0.07987263239/(Foam::pow(moisture, 0.78)))/(Foam::pow(moisture,1.78)) * pSat;



// ************************************************************************* //
