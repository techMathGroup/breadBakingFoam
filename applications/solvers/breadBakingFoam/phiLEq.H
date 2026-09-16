// -- surface liquid flux by Lucas
tmp<volVectorField> jL = -dKoeffLucas * ((Finv.T() & fvc::grad(awPSat)));
surfaceScalarField jLSurf = fvc::interpolate(J * (Finv & jL)) & physics->mesh().Sf();

// -- disable flux through boundaries
forAll(jLSurf.boundaryField(), patchI)
{
    if (jLSurf.boundaryField()[patchI].type() == "calculated")
    {
        scalarField& jLSurfP = jLSurf.boundaryFieldRef()[patchI];
        forAll(jLSurfP, faceI)
        {
            jLSurfP[faceI] = 0;
        }
    }
}

liquidOutFlux = fvc::div(jLSurf);

// -- conservation of the liquid water
fvScalarMatrix concEqL
(
    // -- accumulation
    alphaD0 * rhoSScalar * fvm::ddt(moisture)

    ==

    // -- evaporation source
    - J * evRate

    // -- liquid flux
    - liquidOutFlux
);

concEqL.relax();

phiLResidual = concEqL.solve().initialResidual();

moisture.correctBoundaryConditions();

moisture.relax();



// ************************************************************************* //