// -- saturated water pressure
pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));
pSat.correctBoundaryConditions();


// -- Oswin model for water activity -- legacy
// a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
// a.correctBoundaryConditions();
// aw = Foam::pow(moisture * 100.0 / a, invNCoeff) / (Foam::pow(moisture * 100.0 / a, invNCoeff) + 1.0); 

// -- water activity from Lostie Et al. 2002
aw = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture, 0.78))); 
aw.correctBoundaryConditions();

// -- calculation of free effective surface
// volScalarField alphaL = alphaD * omegaL * rhoD / rhoL;
volScalarField alphaL = alphaD * moisture * rhoSScalar / rhoL;
Sl = alphaL / (alphaL + alphaG);
availSurf = -16 * (Foam::pow(Sl - 0.5, 4)) + 1;
// availSurf =  1;
availSurf.correctBoundaryConditions();
Sl.correctBoundaryConditions();

alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.1415 + 0.5;
// alpha = - Foam::atan((moisture-0.1) / 0.01) / 3.1415 + 0.5;
// alpha = - Foam::atan((Sl-0.04) / 0.001) / 3.1415 + 0.5;
alpha.correctBoundaryConditions();

kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;
kMPC.correctBoundaryConditions();

awPSat = aw * pSat;
awPSat.correctBoundaryConditions();
pV = pG * yV;
pV.correctBoundaryConditions();


// evaporation rate
evRate = kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * (awPSat - pV) * molMV;
evRate.correctBoundaryConditions();

// -- evaporation enthalpy -- legacy
// volScalarField tau = 1 - T / dummyT / 647.126;
// volScalarField dHT = - dimensionedScalar("dummyDH", dH.dimensions(), 1)
//                     * (2419.53 * Foam::pow(tau, 0.333)
//                     - 795.21 * Foam::pow(tau,0.667)
//                     + 6102.03 * tau
//                     - 23771.6 * Foam::pow(tau, 1.667)
//                     + 19241 * Foam::pow(tau, 1.833)) * 1000;

// dimensionedScalar dHConst("dummyDH", dH.dimensions(), -2430492.1);
// dimensionedScalar dHConst("dummyDH", dH.dimensions(), 2257362.5);
// dimensionedScalar dHConst("dummyDH", dH.dimensions(), 2300000);
// dimensionedScalar dHS("dummyDH", dH.dimensions(), 20.7e3);
dimensionedScalar dHS("dummyDH", dH.dimensions(), 573e3);
// ************************************************************************* //