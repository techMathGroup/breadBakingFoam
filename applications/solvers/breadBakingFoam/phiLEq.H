surfaceScalarField SfMag = physics->mesh().magSf();

jLCorr = fvc::interpolate(Finv & jL) & physics->mesh().Sf();
jLCorrIntoLiquid =  0*fvc::interpolate(Finv & jL) & physics->mesh().Sf();
surfaceScalarField omegaVSurf = fvc::interpolate(omegaV);
surfaceScalarField jLCorrOld = jLCorr.prevIter();
surfaceScalarField jLCorrIntoLiquidOld = jLCorrIntoLiquid.prevIter();
surfaceScalarField jPgCorr = fvc::interpolate(jGTilda & fvc::grad(pG)) & physics->mesh().Sf();

// -- zero flux everywhere except BC
// forAll(jPgCorr, faceI)
// {
//     jPgCorr[faceI] = 0;
// }

if (zhang)
{
    forAll(jLCorr.boundaryField(), patchI)
    {
        if (jLCorr.boundaryField()[patchI].type() == "calculated")
        {
            scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
            scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
            scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
            scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
            scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
            scalarField& jLCorrIntoLiquidP = jLCorrIntoLiquid.boundaryFieldRef()[patchI];
            scalarField& jPgCorrP = jPgCorr.boundaryFieldRef()[patchI];
            // scalarField& MgP = Mg.boundaryFieldRef()[patchI];
            // scalarField& TP = T.boundaryFieldRef()[patchI];
            // scalarField& pGP = pG.boundaryFieldRef()[patchI];


            forAll(jLCorrP, faceI)
            {
                dKoeffLucasP[faceI] = 0;
                if (omegaVP[faceI] < SMALL)
                {
                    // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
                    jLCorrIntoLiquidP[faceI] = 0;
                    jLCorrP[faceI] = 0;
                }
                if (omegaVP[faceI] >= SMALL)
                {
                    jLCorrIntoLiquidP[faceI] = 0;
                    // jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
                    jLCorrP[faceI] = 0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
                }
                // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
            }
        }
    }
}
else
{
    forAll(jLCorr.boundaryField(), patchI)
    {
        // if (jLCorr.boundaryField()[patchI].type() == "calculated")
        if (patchI == 0)
        {
            scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
            scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
            scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
            scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
            scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
            scalarField& jLCorrIntoLiquidP = jLCorrIntoLiquid.boundaryFieldRef()[patchI];
            scalarField& jPgCorrP = jPgCorr.boundaryFieldRef()[patchI];

            forAll(jLCorrP, faceI)
            {
                dKoeffLucasP[faceI] = 0;
                if (omegaVP[faceI] < SMALL)
                {
                    // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
                    jLCorrIntoLiquidP[faceI] = 0;
                    jLCorrP[faceI] = 0;
                }
                if (omegaVP[faceI] >= SMALL)
                {
                    jLCorrIntoLiquidP[faceI] = 0;
                    jLCorrP[faceI] = 0.006 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
                }
                // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
            }
        }

        if (patchI == 3)
        {
            scalarField& omegaVP = omegaV.boundaryFieldRef()[patchI];
            scalarField& rhoGP = rhoG.boundaryFieldRef()[patchI];
            scalarField& magSfP = SfMag.boundaryFieldRef()[patchI];
            scalarField& dKoeffLucasP = dKoeffLucas.boundaryFieldRef()[patchI];
            scalarField& jLCorrP = jLCorr.boundaryFieldRef()[patchI];
            scalarField& jLCorrIntoLiquidP = jLCorrIntoLiquid.boundaryFieldRef()[patchI];
            scalarField& jPgCorrP = jPgCorr.boundaryFieldRef()[patchI];


            forAll(jLCorrP, faceI)
            {
                dKoeffLucasP[faceI] = 0;
                if (omegaVP[faceI] < SMALL)
                {
                    // jLCorrIntoLiquidP[faceI] = -0.01 * rhoGP[faceI] * (omegaVP[faceI]) * magSfP[faceI];
                    jLCorrIntoLiquidP[faceI] = 0;
                    jLCorrP[faceI] = 0;
                }
                if (omegaVP[faceI] >= SMALL)
                {
                    jLCorrIntoLiquidP[faceI] = 0;
                    jLCorrP[faceI] = 0.0006 * rhoGP[faceI] * (omegaVP[faceI] - omegaVInfTable_(runTime.time().timeOutputValue())) * magSfP[faceI] - jPgCorrP[faceI] * omegaVP[faceI];
                }
                // jPgCorrP[faceI] = 0.01 * MgP[faceI] / univR.value() * (pGP[faceI] / TP[faceI] - 1e5 / 463) * magSfP[faceI];
            }
        }
    }
}



// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    alphaD0 * fvm::ddt(rhoD)

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


    // - fvc::div(fvc::interpolate(J) * jLCorr)
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

// ************************************************************************* //