surfaceScalarField SfMag = physics->mesh().magSf();

jLCorr = fvc::interpolate(Finv & jL) & physics->mesh().Sf();


forAll(jLCorr.boundaryField(), patchI)
{
    if (jLCorr.boundaryField()[patchI].type() == "calculated")
    {
        scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
        forAll(jLCorrP, faceI)
        {
            jLCorrP[faceI] = 0;
        }
    }
}

// if (zhang)
// {
//     forAll(jLCorr.boundaryField(), patchI)
//     {
//         if (jLCorr.boundaryField()[patchI].type() == "calculated")
//         {
//             scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
//             scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
//             scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
//             scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
//             scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];

//             forAll(jLCorrP, faceI)
//             {
//                 dKoeffLucasP[faceI] = 0;
//                 if (omegaVP[faceI] < SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     jLCorrP[faceI] = 0;
//                 }
//                 if (omegaVP[faceI] >= SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
//                     jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                 }
//                 // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
//             }
//         }
//     }
// }
// else
// {
//     forAll(jLCorr.boundaryField(), patchI)
//     {
//         // if (jLCorr.boundaryField()[patchI].type() == "calculated")
//         if (patchI == 0)
//         {
//             scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
//             scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
//             scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
//             scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
//             scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
//             // scalarField& jLCorrIntoLiquidP = jLCorrIntoLiquid.boundaryFieldRef()[patchI];
//             // scalarField& jPgCorrP = jPgCorr.boundaryFieldRef()[patchI];

//             forAll(jLCorrP, faceI)
//             {
//                 dKoeffLucasP[faceI] = 0;
//                 if (omegaVP[faceI] < SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                     // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
//                     jLCorrP[faceI] = 0;
//                 }
//                 if (omegaVP[faceI] >= SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
//                     jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI];
//                     // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI];
//                 }
//                 // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
//             }
//         }

//         if (patchI == 3)
//         {
//             scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
//             // scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
//             scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
//             scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
//             scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
//             // scalarField& jLCorrIntoLiquidP = jLCorrIntoLiquid.boundaryFieldRef()[patchI];
//             // scalarField& jPgCorrP = jPgCorr.boundaryFieldRef()[patchI];


//             forAll(jLCorrP, faceI)
//             {
//                 dKoeffLucasP[faceI] = 0;
//                 if (omegaVP[faceI] < SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     jLCorrP[faceI] = 0;
//                 }
//                 if (omegaVP[faceI] >= SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     jLCorrP[faceI] = 0;
//                     // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
//                     // jLCorrP[faceI] = - jPgCorrP[faceI] * omegaVP[faceI];
//                 }
//                 // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
//             }
//         }
//     }
// }


// if (zhang)
// {
//     forAll(jLCorr.boundaryField(), patchI)
//     {
//         if (jLCorr.boundaryField()[patchI].type() == "calculated")
//         {
//             scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
//             scalarField& alphaGP = alphaG.boundaryFieldRef()[patchI];
//             scalarField& yVP = yV.boundaryFieldRef()[patchI];
//             scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
//             scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
//             scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
//             scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
//             scalarField& alphaP = alpha.boundaryFieldRef()[patchI];

//             forAll(jLCorrP, faceI)
//             {
//                 // dKoeffLucasP[faceI] = 0;
//                 if (omegaVP[faceI] < SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     jLCorrP[faceI] = 0;
//                 }
//                 if (omegaVP[faceI] >= SMALL)
//                 {
//                     // jLCorrIntoLiquidP[faceI] = 0;
//                     // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
//                     jLCorrP[faceI] = (1-alphaP[faceI]) * 0.01 * rhoGP[faceI] * (omegaVP[faceI] * alphaGP[faceI]) * magSfP[faceI];
//                     // jLCorrP[faceI] = (1-alphaP[faceI]) * 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
//                     // jLCorrP[faceI] = (1-alphaP[faceI]) * 0.01 * rhoGP[faceI] * (yVP[faceI]) * magSfP[faceI];
//                     // jLCorrP[faceI] = (1-alphaP[faceI]) * 0.01 * rhoGP[faceI] * (yVP[faceI]) * magSfP[faceI];
//                 }
//                 // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
//             }
//         }
//     }
// }

odtokMoisture = fvc::div(fvc::interpolate(J) * jLCorr);

// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    // alphaD0 * fvm::ddt(rhoD)
    alphaD0 * rhoSScalar * fvm::ddt(moisture)
    // rhoSScalar * fvm::ddt(alphaD, J, moisture)

    ==

    // - fvc::laplacian(J * Finv & Finv.T() & jGTilda * (1 - alpha) * alphaG * DG / univR / T * molMV, pSat)

    // -- evaporation source
    - J * evRate

    // - J * mCO2

    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc * aw / univR / T * molMV, pSat)
    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc / univR * molMV, pSat * aw / T)
    // + fvc::laplacian(J * Finv * alphaG * (1 - alpha) * Dvc / univR / T * molMV, pV)

    // + fvc::laplacian(J * Finv * dKoeffLucas, awPSat)
    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, pV)


    - fvc::div(fvc::interpolate(J) * jLCorr)
    // - fvc::div(fvc::interpolate(J) * jLCorrIntoLiquid)

    // + fvc::laplacian(J * (Finv & Finv.T()) * dKoeffLucas, awPSat)
    // + fvc::div(jCCorr)
    // + fvc::div(J * dKoeffLucas * (Finv & (Finv.T() & fvc::grad(pV))))
    // + fvc::div(J * dKoeffLucas * (Finv & (Finv.T() & fvc::grad(pV))))
    // + fvc::div(J * alphaG * (1 - alpha) * Dvc / univR / T * molMV * (Finv & (Finv.T() & fvc::grad(pV))))
    // + fvc::div(J * dKoeffLucas * (Finv & (Finv.T() & fvc::grad(awPSat))))

    // - fvm::Sp(citatel / (Foam::pow(T / dummyT, 0.5) * dummyT * Foam::pow(zavorka, 2) * nCoef * alphaL), alphaL)
    // + citatel / (Foam::pow(T / dummyT, 0.5) * dummyT * Foam::pow(zavorka, 2) * nCoef)
);

// SolverPerformance<scalar> solverPerfD;


concEqL.relax();

phiLResidual = concEqL.solve().initialResidual();

rhoD.correctBoundaryConditions();

rhoD.relax();

// mCO2Dis = 0.067275 * Foam::exp(-0.00207 * (T / dummyT - 273)) * fvc::ddt(T) * 1e-3 * rhoD * alphaD * omegaL / dummyT;
dimensionedScalar KH0("KH0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.035);
volScalarField kH = KH0 * Foam::exp(2400.0 * (1.0 / (T / dummyT) - 1.0 / 298.15));
dimensionedScalar kCO2("kCO2", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.2);
volScalarField kCO2VSF = kCO2 * availSurf / 8.314 / Foam::pow(T / dummyT, 0.5) * molMC.value();
// dimensionedScalar kCO2("kCO2", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0);
mCO2Dis = kCO2VSF * alphaG * (1e5 * XCO2InL / kH / molMC.value() * dummyP - yC * pG);
pCSat = 1e5 * XCO2InL / kH / molMC.value() * dummyP;

fvScalarMatrix CO2inLEqn
(
    // alphaD0 * fvm::ddt(rhoD, omegaL, XCO2InL)
    alphaD0 * rhoSScalar * fvm::ddt(moisture, XCO2InL)

    ==

    - fvm::Sp(J * kCO2VSF * (1e5 / kH / molMC.value() * dummyP), XCO2InL)
    // - (kCO2 * (1e5 / kH / molMC.value() * rhoD / dummyRho * omegaL * dummyP) * XCO2InL)
    + J * kCO2VSF * yC * pG

    // - mCO2Dis
);

// Info << "Min (XCO2InL): " << min(XCO2InL).value() << ", max (XCO2InL): " << max(XCO2InL).value() << "." << endl;

CO2inLEqn.relax();

CO2inLEqn.solve().initialResidual();

XCO2InL.correctBoundaryConditions();

XCO2InL.relax();

// ************************************************************************* //