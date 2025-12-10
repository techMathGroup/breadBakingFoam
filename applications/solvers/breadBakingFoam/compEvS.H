// -- partial vapor presure
pV = pG * yV;
// pV = dummyP * 1e5 * yV;

// -- different approaches to evRate calculation
// if (withEq == 1)
// {
    // -- saturated water pressure
    pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));
    // pSat = dummyP * 597 * Foam::exp((T / dummyT - 273) / (T / dummyT - 35.86));

    // Zhang (2025)
    // volScalarField A = 15.64 - 0.1 * (T / dummyT - 273);
    // volScalarField B = 0.38 + 1.69e-3 * (T - 273);

    scalar n = 0.38;
    // volScalarField a = Foam::exp(-0.0071 * T / dummyT + 4.5);
    // volScalarField a = Foam::exp(-0.0056 * T / dummyT + 5.5);
    volScalarField a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
    volScalarField aw = Foam::pow(moisture * 100 / a, (1.0 / n)) / (Foam::pow(moisture * 100/ a, (1.0 / n)) + 1); 
    pVEq = pSat * aw;

    volScalarField alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.14 + 0.5;
    volScalarField kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;

    // evRate = kMPC * (1 - alphaL - alphaS) / univR / T * (aw * pSat - pV) * molMV;
    // evRate = kMPC * (1 - alphaL - alphaS) / Foam::pow(J, 0.333) / univR / T * (aw * pSat - pV) * molMV;
    // evRate = kMPC * (1 - alphaL - alphaS) / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * (aw * pSat - pV) * molMV;
    evRate = kMPC / univR / Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 1), 0.5) / dimensionedScalar("dummyT", dimTemperature, 1) * (aw * pSat - pV) * molMV;

    // -- equalibrium water vapor presure (Zhang)
    // volScalarField a = 4463.083750 * moisture * Foam::exp(0.007 * T / dummyT - 4.61) * (Foam::pow(moisture * Foam::exp(0.007 * T / dummyT - 4.61), 0.82482));
    // pVEq = pSat * (1 / (1 + 10 / a)); // real
    // pVEq = pSat * (1 / (1 + 1 / a)); // in Zhang
    // scalar a = -0.0056;
    // scalar b = 5.5;
    // // scalar a = -0.0056;
    // // scalar b = 5.5;
    // volScalarField aw = 1 / (Foam::pow(moisture / Foam::exp(a * T / dummyT + b), -2.631) + 1);


    // -- Oswin model 
    // scalar n = 0.38;
    // volScalarField a = Foam::exp(-0.0071 * T / dummyT + 4.5);
    // // volScalarField a = Foam::exp(-0.0056 * T / dummyT + 5.5);
    // volScalarField aw = Foam::pow(moisture * 100 / a, (1.0 / n)) / (Foam::pow(moisture * 100/ a, (1.0 / n)) + 1); 
    // pVEq = pSat * aw;

    // Info << "min aw " << min(aw).value() << " max aw " << max(aw).value() << endl;

    // // -- calculation of the evaporation source
    // evRate = intBoil * alphaL / (1-alphaL) * molMC / univR / dimensionedScalar("dummyT", dimTemperature, 300) * (pVEq - pV);
    // evRate = intBoil * (1 - alphaS) * molMC / univR / dimensionedScalar("dummyT", dimTemperature, 300) * (pVEq - pV);
    // evRate = intBoil * (1 - alphaS) * molMC / univR / T * (pVEq - pV);
// }
// -- from (M. Al-Nasser Numerical simulation of bread baking in a convection oven 2021)
// else
// {
//     // -- saturated temperature
//     forAll(TSat, cellI)
//     {
//         if (pV[cellI] > SMALL)
//         {
//             TSat[cellI] = (-2.4432439345097229911e7+2.30650e5*Foam::log(pV[cellI]))/(-1.1598177639491068526e5+5000*Foam::log( pV[cellI]));
//         }
//     }

//     evRate = 0 * evRate;
//     forAll(T, cellI)
//     {
//         if (T[cellI] <= TSat[cellI])
//         {
//             evRate[cellI] +=  c.value() * (1 - alphaS[cellI] - alphaL[cellI]) * rhoG[cellI] * omegaV[cellI] * (T[cellI] - TSat[cellI]) / TSat[cellI];
//         }
//         else
//         {
//             evRate[cellI] += d.value() * alphaL[cellI] * rhoL.value() * (T[cellI] - TSat[cellI]) / TSat[cellI];
//         }
//     }

//     evRate = evRate * Foam::pow(TSat / dummyT, -0.5);
// }

// ************************************************************************* //