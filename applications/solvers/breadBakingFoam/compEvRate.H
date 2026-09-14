// -- saturated water pressure
pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));
pSat.correctBoundaryConditions();


// -- Oswin model for water activity -- legacy
a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
a.correctBoundaryConditions();
aw = Foam::pow(moisture * 100.0 / a, invNCoeff) / (Foam::pow(moisture * 100.0 / a, invNCoeff) + 1.0); 

// Info << "moisture min: " << min(moisture).value() << ", max: " << max(moisture).value() << endl;
volScalarField alphaL = alphaD * moisture * rhoSScalar / rhoL;
Sl = alphaL / (alphaL + alphaG);



forAll(moisture, i)
{
    if (moisture[i] < SMALL)
    {
        // aw[i] = 0.0;
        availSurf[i] = 0.0;
    }
    else
    {
        // aw[i] = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture[i], 0.78))); 
        availSurf[i] = -16 * (Foam::pow(Sl[i] - 0.5, 4)) + 1;
    }
}

forAll(aw.boundaryField(), patchI)
{
    scalarField& awP = aw.boundaryFieldRef()[patchI];
    scalarField& availSurfP = availSurf.boundaryFieldRef()[patchI];
    scalarField& moistureP = moisture.boundaryFieldRef()[patchI];
    forAll(awP, faceI)
    {
        if (moistureP[faceI] < SMALL)
        {
            // awP[faceI] = 0.0;
            availSurfP[faceI] = 0.0;
        }
        else
        {
            // awP[faceI] = Foam::exp(-2.9 / (Foam::pow(100.0 * moistureP[faceI], 0.78))); 
            availSurfP[faceI] = -16 * (Foam::pow(Sl[faceI] - 0.5, 4)) + 1;
            // availSurfP[faceI] = 1 - Foam::mag(2 * (Sl[faceI] - 0.5));
        }
    }
}

// -- water activity from Lostie Et al. 2002
// aw = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture, 0.78))); 
// aw.correctBoundaryConditions();

// -- calculation of free effective surface
// volScalarField alphaL = alphaD * omegaL * rhoD / rhoL;
// volScalarField alphaL = alphaD * moisture * rhoSScalar / rhoL;
// Sl = alphaL / (alphaL + alphaG);
// availSurf = -16 * (Foam::pow(Sl - 0.5, 4)) + 1;
availSurf =  1;
availSurf.correctBoundaryConditions();
// Sl.correctBoundaryConditions();

alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.1415 + 0.5;
// alpha = - Foam::atan((moisture-0.1) / 0.005) / 3.1415 + 0.5;
// alpha = - Foam::atan((moisture-moistureMin) / 0.001) / 3.1415 + 0.5;
// alpha = 1.0 - Foam::exp(- 10.0 * (0.5 - moisture));
// forAll(alpha, i)
// {
//     if (moisture[i] > 0.5)
//     {
//         // aw[i] = 0.0;
//         alpha[i] = SMALL;
//     }
//     else
//     {
//         // aw[i] = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture[i], 0.78))); 
//         // availSurf[i] = -16 * (Foam::pow(Sl[i] - 0.5, 4)) + 1;
//         alpha[i] = 1.0 - Foam::exp(- 10.0 * (0.5 - moisture[i]));
//     }
// }

forAll(alpha.boundaryField(), patchI)
{
    if (jLCorr.boundaryField()[patchI].type() == "calculated")
    {
        scalarField& alphaP = alpha.boundaryFieldRef()[patchI];
        // scalarField& moistureP = availSurf.boundaryFieldRef()[patchI];
        scalarField& moistureP = moisture.boundaryFieldRef()[patchI];
        forAll(alphaP, faceI)
        {
            // if (moistureP[faceI]  > 0.5)
            // {
            //     // awP[faceI] = 0.0;
            //     alphaP[faceI] = SMALL;
            // }
            // else
            // {
            //     // aw[i] = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture[i], 0.78))); 
            //     // availSurf[i] = -16 * (Foam::pow(Sl[i] - 0.5, 4)) + 1;
            //     // alphaP[faceI] = 1.0 - Foam::exp(- 10.0 * (0.5 - moistureP[faceI]));
            //     alphaP[faceI] = - Foam::atan((moistureP[faceI]-moistureMin) / 0.001) / 3.1415 + 0.5;
            // }
            alphaP[faceI] = 1 - SMALL;
        }
    }
}
// alpha = - Foam::atan((Sl-0.035) / 0.001) / 3.1415 + 0.5;
// alpha.correctBoundaryConditions();

kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;
kMPC.correctBoundaryConditions();

awPSat = aw * pSat;
awPSat.correctBoundaryConditions();
pV = pG * yV;
pV.correctBoundaryConditions();

// volScalarField availMoist = Foam::atan((moisture-moistureMin) / 0.001) / 3.1415 + 0.5;


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