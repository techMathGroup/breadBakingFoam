// -- saturated water pressure
pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));

// -- vapor molar fraction
yV = omegaV / molMV * Mg;



// -- Oswin model for water activity
volScalarField a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
// volScalarField availMois = max(moisture - 0.3, SMALL);
volScalarField aw = Foam::pow(moisture * 100 / a, (1.0 / nCoef)) / (Foam::pow(moisture * 100/ a, (1.0 / nCoef)) + 1); 
// volScalarField aw = Foam::pow(availMois * 100 / a, (1.0 / nCoef)) / (Foam::pow(availMois * 100/ a, (1.0 / nCoef)) + 1); 

alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.14 + 0.5;
volScalarField kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;

awPSat = aw * pSat;
pV = pG * yV;

// dimensionedScalar Slr("Slr", dimless, 0.05);
// dimensionedScalar Sgr("Sgr", dimless, 0.02);
// volScalarField Sl = alphaL / (alphaL + alphaG);
// volScalarField Se = (Sl - Slr) / (1 - Slr - Sgr);
// availSurf = 4 * Se * (1 - Se);

// evaporation rate
evRate = kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * (aw * pSat - pG * yV) * molMV;

// -- evaporation enthalpy
volScalarField tau = 1 - T / dummyT / 647.126;
volScalarField dHT = - dimensionedScalar("dummyDH", dH.dimensions(), 1)
                    * (2419.53 * Foam::pow(tau, 0.333)
                    - 795.21 * Foam::pow(tau,0.667)
                    + 6102.03 * tau
                    - 23771.6 * Foam::pow(tau, 1.667)
                    + 19241 * Foam::pow(tau, 1.833)) * 1000;

// dimensionedScalar dHConst("dummyDH", dH.dimensions(), -2430492.1);
// dimensionedScalar dHConst("dummyDH", dH.dimensions(), 2257362.5);
dimensionedScalar dHConst("dummyDH", dH.dimensions(), 2300000);
// dimensionedScalar dHS("dummyDH", dH.dimensions(), 20.7e3);
dimensionedScalar dHS("dummyDH", dH.dimensions(), 573e3);


// jDiffNow = - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG)) + (jDVTilda & fvc::grad(Mg) / Mg)) + (jDVTilda & fvc::grad(omegaV))));
// jDiffNow.correctBoundaryConditions();
// delta = 1 - Foam::pow((alphaG), 1.0/3.0);   
// delta.correctBoundaryConditions();
// pV.correctBoundaryConditions();

// // volScalarField dPSatdT = 3.81644e6*Foam::exp((81.443*T/dummyT - 22839.16559)/(5.*T/dummyT - 230.65))/(Foam::pow(T/dummyT - 46.13, 2)) * dummyP / dummyT;

// jDiffNew = - J * (Finv & (omegaV * ((jGTilda & fvc::grad(pG))) + ( (alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / DG) + delta * dHConst / lambdaS * rhoG * univR * omegaV / molMV * aw )) * Finv.T() & fvc::grad(pV))));
// jDiffNew.correctBoundaryConditions();
// ************************************************************************* //