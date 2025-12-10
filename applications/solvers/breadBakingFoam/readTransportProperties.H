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
scalar withDeformation;
transportProperties.subDict("genProps").readEntry("univR", univR); // -- universal gas constant
transportProperties.readEntry("withDeformation", withDeformation);
// transportProperties.readEntry("withEq", withEq);
transportProperties.readEntry("permGLViscG",permGLViscGSc);

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

// dimensionedScalar kappaL (dimensionedScalar("kappaL", kappaA.dimensions(), 0.643));
// dimensionedScalar kappaL (dimensionedScalar("kappaL", kappaA.dimensions(), 1.5));
// dimensionedScalar kappaG (dimensionedScalar("kappaG", kappaA.dimensions(), 0.0189));
// dimensionedScalar CpL("CpL", dimEnergy/dimMass/dimTemperature, 4200);
// dimensionedScalar CpVapor ("CpVapor", CpL.dimensions(), 1930.55);

// -- reaction
dimensionedScalar dH, R0, kMPCOpen, kMPCClosed, TOpen, deltaTOpen, Tm, deltaT, evCoef1, evCoef2;
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("dH", dH);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("TOpen", TOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("deltaTOpen", deltaTOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("kMPCOpen", kMPCOpen);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("kMPCClosed", kMPCClosed);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("evCoef1", evCoef1);
reactiveProperties.subDict("reactProps").subDict("evaporation").readEntry("evCoef2", evCoef2);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("Tm", Tm);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("R0", R0);
reactiveProperties.subDict("reactProps").subDict("fermentation").readEntry("deltaT", deltaT);
// dimensionedScalar Tm ("Tm", dimTemperature, 313);
// dimensionedScalar deltaT ("Tm", dimTemperature, 10);

// scalar TOpen = 75;
// scalar deltaTOpen = 10;
// dimensionedScalar kMPCOpen ("kMPCOpen", dimless / dimTime, 0.5);
// dimensionedScalar kMPCClosed ("kMPCClosed", dimless / dimTime, 10);

// -- dummies
dimensionedScalar dummyT (dimensionedScalar("dummyT", dimTemperature, 1));
dimensionedScalar dummyP (dimensionedScalar("dummyp", dimPressure, 1));
dimensionedScalar dummyCp (dimensionedScalar("dummyp", CpS.dimensions(), 1));
dimensionedScalar dummyLamdba (dimensionedScalar("dummyLamdba", lambdaS.dimensions(), 1));

double pResidual = 1;
double omegaVResidual = 1;
double phiLResidual = 1;
double TResidual = 1;

// ************************************************************************* //