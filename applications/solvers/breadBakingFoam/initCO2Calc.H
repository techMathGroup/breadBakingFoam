// -- CO2 bookkeeping (initialized on first timestep)
static bool co2Init = false;
static double initialBreadMass = 0.0;  // total mass (solid + liquid + gas) at t0
static double cumCO2 = 0.0;            // cumulative CO2 produced [kg]

if (!co2Init)
{
    const dimensionedScalar breadMass =
        fvc::domainIntegrate
        (
            J *
            (
                alphaD * rhoD * omegaS
                + alphaD * rhoD * omegaL
                + alphaG * rhoG
            )
        );

    initialBreadMass = breadMass.value();
    co2Init = true;
}

const double dtVal = runTime.deltaT().value();
