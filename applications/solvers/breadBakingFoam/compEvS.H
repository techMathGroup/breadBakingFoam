// -- saturated water pressure
pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));

// -- vapor molar fraction
yV = omegaV / molMV * Mg;

// -- Oswin model for water activity
volScalarField a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
volScalarField aw = Foam::pow(moisture * 100 / a, (1.0 / nCoef)) / (Foam::pow(moisture * 100/ a, (1.0 / nCoef)) + 1); 

alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.14 + 0.5;
volScalarField kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;

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
dimensionedScalar dHConst("dummyDH", dH.dimensions(), 2257362.5);
dimensionedScalar dHS("dummyDH", dH.dimensions(), 20.7e3);

// ************************************************************************* //