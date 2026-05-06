Info<< "Reading transport and thermophysical properties\n" << endl; 

// -- transportProperties dictionary
IOdictionary transportProperties
(
    IOobject
    (
	"transportProperties",    // dictionary name
	runTime.constant(),     // dict is found in "constant"
	physics->mesh(),                   // registry for the dict
	IOobject::MUST_READ,    // must exist, otherwise failure
	IOobject::NO_WRITE      // dict is only read by the solver
    )
);

// -- thermophysicalProperties dictionary
IOdictionary thermophysicalProperties
(
    IOobject
    (
	"thermophysicalProperties",    // dictionary name
	runTime.constant(),     // dict is found in "constant"
	physics->mesh(),                   // registry for the dict
	IOobject::MUST_READ,    // must exist, otherwise failure
	IOobject::NO_WRITE      // dict is only read by the solver
    )
);

// -- reactiveProperties dictionary
IOdictionary reactiveProperties
(
    IOobject
    (
        "reactiveProperties",      // dictionary name
        runTime.constant(),            // dict is found in "constant"
        physics->mesh(),                      // registry for the dict
        IOobject::MUST_READ,        // must exist, otherwise failure
        IOobject::NO_WRITE          // dict is only read by the solver
    )
);

// -- model transport parameters 
dimensionedScalar univR, permGLViscGSc, kappaS;
scalar withDeformation, tort;
transportProperties.subDict("genProps").readEntry("univR", univR); // -- universal gas constant
transportProperties.readEntry("withDeformation", withDeformation);
transportProperties.readEntry("permGLViscG",permGLViscGSc);
transportProperties.readEntry("tort",tort);

// -- model thermophysical parameters 
dimensionedScalar rhoS, rhoL, CpS, CpL, CpVapor, CpCO2, molMC, molMV, lambdaS, lambdaL, DG;
thermophysicalProperties.subDict("solid").readEntry("rho", rhoS);
thermophysicalProperties.subDict("liquid").readEntry("rho", rhoL);
thermophysicalProperties.subDict("solid").readEntry("Cp", CpS);
thermophysicalProperties.subDict("liquid").readEntry("Cp", CpL);
thermophysicalProperties.subDict("vapor").readEntry("Cp", CpVapor);
thermophysicalProperties.subDict("CO2").readEntry("Cp", CpCO2);
thermophysicalProperties.subDict("solid").readEntry("lambda", lambdaS);
thermophysicalProperties.subDict("liquid").readEntry("lambda", lambdaL);
thermophysicalProperties.subDict("vapor").readEntry("molM", molMV);
thermophysicalProperties.subDict("CO2").readEntry("molM", molMC);
thermophysicalProperties.subDict("mixture").subDict("transport").readEntry("D", DG);

dimensionedScalar molMAir (dimensionedScalar("molMAir", dimMass/dimMoles, 0.02897));
dimensionedScalar CpAir = dimensionedScalar("CpAir", dimEnergy/dimMass/dimTemperature, 1005);

// -- reaction
dimensionedScalar dH, R0, kMPCOpen, kMPCClosed, TOpen, deltaTOpen, Tm, deltaT, evCoef1, evCoef2;
scalar nCoef, TGelatEv, alphaKept;
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("dH", dH);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("TOpen", TOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("deltaTOpen", deltaTOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("kMPCOpen", kMPCOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("kMPCClosed", kMPCClosed);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("evCoef1", evCoef1);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("evCoef2", evCoef2);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("nCoef", nCoef);
reactiveProperties.subDict("reactProps").subDict("gelatization").readEntry("TGelEv", TGelatEv);
reactiveProperties.subDict("reactProps").subDict("gelatization").readEntry("alphaKept", alphaKept);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("Tm", Tm);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("R0", R0);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("deltaT", deltaT);

// -- dummies
dimensionedScalar dummyT (dimensionedScalar("dummyT", dimTemperature, 1));
dimensionedScalar dummyP (dimensionedScalar("dummyp", dimPressure, 1));
dimensionedScalar dummyCp (dimensionedScalar("dummyp", CpS.dimensions(), 1));
dimensionedScalar dummyLamdba (dimensionedScalar("dummyLamdba", lambdaS.dimensions(), 1));
dimensionedScalar dummyRho (dimensionedScalar("dummyRho", rhoS.dimensions(), 1));

double pResidual = 1;
double omegaVResidual = 1;
double omegaCResidual = 1;
double phiLResidual = 1;
double TResidual = 1;

// ************************************************************************* //