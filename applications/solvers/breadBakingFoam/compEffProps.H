// -- calculation of the Darcy permeability, effective difusivity and effective heat conductivity

Dvc = dimensionedScalar("Dvc", DG.dimensions(), 2.15e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dva = dimensionedScalar("Dva", DG.dimensions(), 2.65e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dca = dimensionedScalar("Dca", DG.dimensions(), 1.67e-5) * Foam::pow(T / dimensionedScalar("dummyT", dimTemperature, 300), 1.75); // -- mass difusivity at correct temperature
Dvc.correctBoundaryConditions();
Dva.correctBoundaryConditions();
Dca.correctBoundaryConditions();

nuG = 1.52451e-6 * Foam::pow(T / dummyT, 1.5) / (T / dummyT + 232.113); // -- gas viscosity at correct temperature
nuG.correctBoundaryConditions();

tortField = tortOpen * alpha + (1 - alpha) * tortClosed;
tortField.correctBoundaryConditions();

// -- effective diffusivity
// DEff = D * alphaG / tort;
// DEff = D * alphaG / tortField * (Foam::atan((moisture-0.35)/0.15)/3.14 + 0.5);

DEffvM = alphaG * alpha / tortOpen * 1.0 / ((yC / Dvc + yA / Dva) + yV / (1 - omegaV) * (omegaC / Dvc + omegaAir / Dva));
DEffcM = alphaG * alpha / tortOpen * 1.0 / ((yV / Dvc + yA / Dca) + yC / (1 - omegaC) * (omegaV / Dvc + omegaAir / Dca));
DEffaM = alphaG * alpha / tortOpen * 1.0 / ((yC / Dca + yV / Dva) + yA / (1 - omegaAir) * (omegaC / Dca + omegaV / Dva));
DEffvM.correctBoundaryConditions();
DEffcM.correctBoundaryConditions();
DEffaM.correctBoundaryConditions();

// -- permeability computation
dimensionedScalar kBaked("k_baked", permGLViscGSc.dimensions(), 1e-11);   

// -- logarithmic interpolation constants
scalar logKRaw = Foam::log(permGLViscGSc.value());
scalar logKBaked = Foam::log(kBaked.value());

// -- permiability MPC Zhang 2025
// permGLViscG = alpha * permGLViscGSc / nuG * Foam::pow(alphaG, 1.34); 
permGLViscG =  Foam::exp(logKRaw + alpha * (logKBaked - logKRaw)) * dimensionedScalar("permGLViscGSc", permGLViscGSc.dimensions(), 1.0) / nuG; 
// permGLViscG =  permGLViscGSc / nuG; 
// permGLViscG =  permGLViscGSc / nuG * Foam::pow(alphaG, 1.34); 
// permGLViscG =  alpha * permGLViscGSc / nuG ; 
permGLViscG.correctBoundaryConditions();

// -- effective heat conductivity https://www.sciencedirect.com/science/article/pii/S0960308508000333
// lambdaEff = dummyLamdba * (0.421 - 0.134 * (1 - alphaS - alphaL) - 0.219 * Foam::pow(1 - alphaS - alphaL, 0.5));
// lambdaEff = dummyLamdba * (lambdaS.value() - 0.134 * (1 - alphaS - alphaL) - 0.219 * Foam::pow(1 - alphaS - alphaL, 0.5));
// lambdaEff = dummyLamdba * (Foam::atan((moisture-0.35)/0.15)/3.14 + 0.5)  * (lambdaS.value() - 0.134 * alphaG - 0.219 * Foam::pow(alphaG, 0.5));
// lambdaEff = dummyLamdba * (Foam::atan((moisture-0.35)/0.15)/3.14 + 0.5)  * (lambdaS.value());
lambdaEff = dummyLamdba *  (lambdaS.value() - 0.134 * alphaG - 0.219 * Foam::pow(alphaG, 0.5));

// -- purlis lambda alternative
// forAll(lambdaEff, cellI)
// {
//     if (T[cellI] > 373)
//     {
//         lambdaEff[cellI] = 0.2;
//     }
//     else
//     {
//         lambdaEff[cellI] = 0.9 / (1.0 + Foam::exp(-0.1 * (T[cellI] - 353.16))) + 0.2;
//     }
// }
// lambdaEff.correctBoundaryConditions();
// lambdaEff = dummyLamdba *  (lambdaS.value());
// lambdaEff = dummyLamdba *  (lambdaS.value() * alphaD + 0.02 * alphaG);

// -- lucas lambda alternative
// volScalarField phi = (moisture - 0.11111) / (0.8518 - 0.11111);
// lambdaEff = dummyLamdba *  (0.423 * phi + 0.18 * (1-phi));
// lambdaEff = dummyLamdba *  (0.37 * phi + 0.18 * (1-phi));
// lambdaEff = dummyLamdba *  0.25;

lambdaEff.correctBoundaryConditions();

// -- calculation of variables for Lucas closed-cell vapor flux
dPSatdT = 3.81644e6*Foam::exp((81.443*T/dummyT - 22839.16559)/(5.*T/dummyT - 230.65))/(Foam::pow(T/dummyT - 46.13, 2)) * dummyP / dummyT;
dPSatdT.correctBoundaryConditions();

volScalarField dPVdT = rhoG * univR / Mg * yV;
dPVdT.correctBoundaryConditions();

forAll(dPWdX, cellI)
{
    if (moisture[cellI] > SMALL)
    {
        dPWdX[cellI] = ((Foam::pow(100,(2.0/nCoef)))*(Foam::pow((moisture[cellI]/a[cellI]),(2.0/nCoef))) - Foam::pow(10000,(1.0/nCoef))*Foam::pow((moisture[cellI]/a[cellI]),(2.0/nCoef)) + Foam::pow(100,(1.0/nCoef))*Foam::pow((moisture[cellI]/a[cellI]),(1.0/nCoef))/(nCoef*moisture[cellI]*Foam::pow((Foam::pow(100,(1.0/nCoef))*Foam::pow((moisture[cellI]/a[cellI]),(1.0/nCoef)) + 1),2))) * pSat[cellI];
    }
    else
    {
        dPWdX[cellI] = 0;
    }
}

forAll(dPWdX.boundaryField(), patchI)
{
    scalarField& dPWdXP = dPWdX.boundaryFieldRef()[patchI];
    scalarField& moistureP = moisture.boundaryFieldRef()[patchI];
    scalarField& aP = a.boundaryFieldRef()[patchI];
    scalarField& pSatP = pSat.boundaryFieldRef()[patchI];

    forAll(dPWdXP, faceI)
    {
        if (moistureP[faceI] > SMALL)
        {
            dPWdXP[faceI] = ((Foam::pow(100,(2.0/nCoef)))*(Foam::pow((moistureP[faceI]/aP[faceI]),(2.0/nCoef))) - Foam::pow(10000,(1.0/nCoef))*Foam::pow((moistureP[faceI]/aP[faceI]),(2.0/nCoef)) + Foam::pow(100,(1.0/nCoef))*Foam::pow((moistureP[faceI]/aP[faceI]),(1.0/nCoef))/(nCoef*moistureP[faceI]*Foam::pow((Foam::pow(100,(1.0/nCoef))*Foam::pow((moistureP[faceI]/aP[faceI]),(1.0/nCoef)) + 1),2))) * pSatP[faceI];
        }
        else
        {
            dPWdXP[faceI] = 0;
        }
    }
}

// -- ratio of cell wall to all thickness
// delta = 1 - Foam::pow((alphaG), 1.0/3.0);
// delta = 1.0;
delta = 0.1;
delta.correctBoundaryConditions();
dKoeffLucas = (1 - alpha) * alphaG / ((1 - delta) * (univR * T / molMV * (1 - yV) / Dvc) + delta * (1.0 / (rhoL * Dl) * dPWdX + dHConst / lambdaS * dPSatdT) );
dKoeffLucas.correctBoundaryConditions();


// volScalarField a = Foam::exp(evCoef1 * T / dummyT + evCoef2);
// volScalarField aw = Foam::pow(moisture * 100 / a, (1.0 / nCoef)) / (Foam::pow(moisture * 100/ a, (1.0 / nCoef)) + 1); 


// Info << "min(dPSatdT) " << min(dPSatdT).value() << endl;
// Info << "max(dPSatdT) " << max(dPSatdT).value() << endl;

// Info << "min(aw) " << (1 - alpha) * alphaG * DEff * aw * (-dHT) / univR / T * dPSatdT * dummyP / dummyT * molMV << endl;
// -- zhang lambda alternative
// volScalarField addToLam = (1 - alpha) * alphaG * D * aw * (-dHT) / univR / T * dPSatdT * dummyP / dummyT * molMV * 5;
// volScalarField addToLam = (1 - alpha) * alphaG * Dvc * aw * (-dHT) / univR / T * dPSatdT * molMV;
// volScalarField addToLam = (1 - alpha) * alphaG * Dvc * aw * (-dHT) / univR / T * dPSatdT * molMV;
// volScalarField addToLam = (1 - alpha) * alphaG * Dvc * (-dHT) / univR / T * dPVdT * molMV;
// volScalarField addToLam = dKoeffLucas * (dHConst) * aw * dPSatdT;

// lambdaEff += addToLam;
// lambdaEff.correctBoundaryConditions();

// Info << "min(addToLam) " << min(addToLam).value() << endl;
// Info << "max(addToLam) " << max(addToLam).value() << endl;
 
// forAll(lambdaEff, cellI)
// {
//     scalar dTepl = 4e-3;
//     scalar yC = 1.0417e-2;
//     vector cellReal = physics->mesh().C()[cellI] + physics->D()[cellI];
//     if (cellReal[0] > 4.45e-2 && cellReal[1] > yC - dTepl/2 && cellReal[1] < yC + dTepl/2)
//     {
//         lambdaEff[cellI] = 0.2;
//     }
// }
// lambdaEff.correctBoundaryConditions();




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