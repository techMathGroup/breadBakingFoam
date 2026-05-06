// -- calculation of the Darcy permeability, effective difusivity and effective heat conductivity

volScalarField D = DG * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
D.correctBoundaryConditions();
nuG = 1.52451e-6 * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113); // -- gas viscosity at correct temperature
nuG.correctBoundaryConditions();

// -- effective diffusivity
DEff = D * alphaG / tort;
DEff.correctBoundaryConditions();

// -- permiability MPC Zhang 2025
permGLViscG = alpha * permGLViscGSc / nuG * Foam::pow(alphaG, 1.34); 
permGLViscG.correctBoundaryConditions();

// -- effective heat conductivity https://www.sciencedirect.com/science/article/pii/S0960308508000333
// lambdaEff = dummyLamdba * (0.421 - 0.134 * (1 - alphaS - alphaL) - 0.219 * Foam::pow(1 - alphaS - alphaL, 0.5));
lambdaEff = dummyLamdba * (lambdaS.value() - 0.134 * (1 - alphaS - alphaL) - 0.219 * Foam::pow(1 - alphaS - alphaL, 0.5));
lambdaEff.correctBoundaryConditions();


// -- legacy from Zhang 2005
// kGForPerm = (1 - 1.1 * (alphaL / (1 - alphaS))); // -- relative permeability (Zhang Transport Processes and Large Deformation During Baking of Bread 2005)
// kGForPerm.correctBoundaryConditions();
// forAll(permGLViscG, cellI)
// {
//     // -- crop for low values of alphaG
//     // if(kGForPerm[cellI] > 1e-4)
//     // if(alphaCrust[cellI] > 0.15)
//     if ((alphaL[cellI] / (1 - alphaS[cellI])) < 0.9)
//     {
//         permGLViscG[cellI] = permGLViscGSc.value() / nuG[cellI] * kGForPerm[cellI];
//     }
//     else
//     {
//         permGLViscG[cellI] = 0;
//     }
//     // {
//         // permGLViscG[cellI] = permGLViscGSc.value() * kGForPerm[cellI] * (1 - alphaS[cellI] - alphaL[cellI]);
//         // DEff[cellI] = D[cellI] * Foam::pow(kGForPerm[cellI] * (1 - alphaS[cellI]), 1.333) * (1 - alphaS[cellI] - alphaL[cellI]);
//         // permGLViscG[cellI] = permGLViscGSc.value() / nuG[cellI] * kGForPerm[cellI];
//         // DEff[cellI] = D[cellI] * Foam::pow(kGForPerm[cellI] * (1 - alphaS[cellI]), 1.333);
//         // DEff[cellI] = D[cellI] * 0.1;
//     // }
//     // else
//     // {
//     //     permGLViscG[cellI] = 1e-11 ;
//     //     // DEff[cellI] = 1e-9;
//     // }
//     // if((1 - alphaL[cellI] - alphaS[cellI]) > 0.75)
//     // {
//     //     permGLViscG[cellI] = permGLViscG[cellI] * 10;
//     // }
// }

// DEff = D * (1 - (Foam::atan(100*(alphaCrust-0.1))/(3.14/2) + 1)/2);

// permGLViscG = permGLViscGSc / nuG * kGForPerm;
// DEff = D * (1 - alphaL - alphaS) * Foam::pow(kGForPerm * (1 - alphaS), 1.333);


// DEff = D * (1 - alphaL - alphaS) / 3.5;
// Info << "min(nuG) " << min(nuG).value() << endl;
// Info << "min(permGLViscG) " << min(permGLViscG).value() << endl;

// permGLViscG = permGLViscGSc / nuG * (Foam::atan(T / dummyT - 323) / 3.14 + 0.5);

// DEff = alphaCrust * dimensionedScalar("smallDEff", DEff.dimensions(), 1e-11) + (1 - alphaCrust) * D;

// permGLViscG = permGLViscGSc / nuG;

// permGLViscG *= (1 - alphaL - alphaS); 

// -- effective heat conductivity 

// -- gas heat conductivity
// volScalarField kappa (kappaA + kappaB * T + kappaC * Foam::pow(T, 2) + kappaD * Foam::pow(T, 3));
// volScalarField lambdaG = dummyLamdba * 2.56008e-3 * Foam::pow(T / dummyT, 0.5 / (1 + 207.009 / T * dummyT));

// -- effective thermal conductivity (different approaches)
// lambdaEff = dummyLamdba * 
//           Foam::pow(lambdaS / dummyLamdba, alphaS) * 
//           Foam::pow(lambdaL / dummyLamdba, alphaL) * 
//           Foam::pow(lambdaG / dummyLamdba, 1 - alphaS - alphaL);

// lambdaEff = (alphaS - alphaS + 1) * lambdaS;
// lambdaEff.correctBoundaryConditions();

// kappaPM = (kappaS * alphaS) + 
//           (kappaL * alphaF * phiL) + 
//           (kappa * alphaF * phiG);
// lambdaEff = (lambdaS * alphaS) + 
//             (lambdaL * alphaL) + 
//             (lambdaG * (1 - alphaS - alphaL));
// lambdaEff = dummyLamdba * (0.03-(0.02/(moisture-0.6)));
// forAll(lambdaEff, cellI)
// {
//     if (lambdaEff[cellI] > 0.3)
//     {
//         lambdaEff[cellI] = 0.3;
//     }
// }



// ************************************************************************* //