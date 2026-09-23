// -- saturated water pressure
pSat = dummyP * 1e3 * Foam::exp(16.2886 - 3816.44 / (T / dummyT - 46.13));
pSat.correctBoundaryConditions();

// -- Oswin model for water activity
moistureSafe = max(moisture, moistureFloor);
a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
a.correctBoundaryConditions();

// -- update water activity
forAll(aw, cellI)
{
    if (moistureSafe[cellI] > SMALL)
    {
        aw[cellI] = Foam::pow(moistureSafe[cellI] * 100.0 / a[cellI], invNCoeff) / (Foam::pow(moistureSafe[cellI] * 100.0 / a[cellI], invNCoeff) + 1.0);    
    } 
    else
    {
        aw[cellI] = SMALL;
    }
}
forAll(aw.boundaryField(), patchI)
{
    scalarField& awP = aw.boundaryFieldRef()[patchI];
    scalarField& moistureP = moisture.boundaryFieldRef()[patchI];
    scalarField& aP = a.boundaryFieldRef()[patchI];
    forAll(awP, faceI)
    {
        if (moistureP[faceI] > SMALL)
        {
            awP[faceI] = Foam::pow(moistureP[faceI] * 100.0 / aP[faceI], invNCoeff) / (Foam::pow(moistureP[faceI] * 100.0 / aP[faceI], invNCoeff) + 1.0); 
        }
        else
        {
            awP[faceI] = SMALL;
        }
    }
}


// -- regime transition function
// alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.1415 + 0.5;
alpha = Foam::atan((T / dummyT - TOpen - 273.15) / (deltaTOpen * 0.125)) / 3.1415 + 0.5;


// -- boundary conditions are opened
forAll(alpha.boundaryField(), patchI)
{
    if (jLSurf.boundaryField()[patchI].type() == "calculated")
    {
        scalarField& alphaP = alpha.boundaryFieldRef()[patchI];
        forAll(alphaP, faceI)
        {
            alphaP[faceI] = 1 - SMALL;
        }
    }
}

// -- evaporation coefficient
kMPC = alpha * kMPCOpen + (1 - alpha) * kMPCClosed;
kMPC.correctBoundaryConditions();

// -- saturated vapor pressure
awPSat = aw * pSat;
awPSat.correctBoundaryConditions();
pV = pG * yV;
pV.correctBoundaryConditions();

// -- evaporation rate
evRate = kMPC * availSurf / univR / Foam::pow(T / dummyT, 0.5) / dummyT * (awPSat - pV) * molMV;
evRate.correctBoundaryConditions();

// -- water activity from Lostie Et al. 2002 (legacyy)
// aw = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture, 0.78))); 
// aw.correctBoundaryConditions();

// -- computation of the available surface -- legacy
// volScalarField alphaL = alphaD * moisture * rhoSScalar / rhoL;
// Sl = alphaL / (alphaL + alphaG);
// forAll(moisture, i)
// {
//     if (moisture[i] < SMALL)
//     {
//         // aw[i] = 0.0;
//         availSurf[i] = 0.0;
//     }
//     else
//     {
//         // aw[i] = Foam::exp(-2.9 / (Foam::pow(100.0 * moisture[i], 0.78))); 
//         availSurf[i] = -16 * (Foam::pow(Sl[i] - 0.5, 4)) + 1;
//     }
// }

// forAll(aw.boundaryField(), patchI)
// {
//     scalarField& awP = aw.boundaryFieldRef()[patchI];
//     scalarField& availSurfP = availSurf.boundaryFieldRef()[patchI];
//     scalarField& moistureP = moisture.boundaryFieldRef()[patchI];
//     forAll(awP, faceI)
//     {
//         if (moistureP[faceI] < SMALL)
//         {
//             // awP[faceI] = 0.0;
//             availSurfP[faceI] = 0.0;
//         }
//         else
//         {
//             // awP[faceI] = Foam::exp(-2.9 / (Foam::pow(100.0 * moistureP[faceI], 0.78))); 
//             availSurfP[faceI] = -16 * (Foam::pow(Sl[faceI] - 0.5, 4)) + 1;
//             // availSurfP[faceI] = 1 - Foam::mag(2 * (Sl[faceI] - 0.5));
//         }
//     }
// }

// ************************************************************************* //