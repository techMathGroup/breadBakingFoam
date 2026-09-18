# 2D our custom experiment
## Case description and setup
This tutorial shows a two-dimensional internal simulation of the bread free proofing and baking in our laboratory oven. External transport is resolved by custom mixed boundary conditions. The tutorial is located in `tutorials/breadAx2D` and can be:
1. run directly as prepared by `Allrun` script in `tutorials/breadAx2DOurExp` folder, or
2. modified and run by `pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py` control script.

The description of the solved equations and variables is in greater detail discussed in https://doi.org/10.14311/TPFM.2025.015. Furthermore in the solver, the solved variables are noted as: 
* `moisture` - relative mass fraction of liquid water with respect to solid mass
* `T` - temperature,
* `pG` - pressure of the gas phase,
* `omegaV` - mass fraction of the water vapors in the gas,
* `omegaC` - mass fraction of the CO2 in the gas, and
* `D` - deformation vector. 

### Geometry and computational mesh description
Geometry is based on our custom experiments conducted in the laboratories at University of Chemistry and Technology. The two breads are simultaniously placed into oven to measure both the temperature evolution at several different places in the bread and the bread weight. 

<img alt="tutBreadCustomExpDescr" src="tutBread2DOurGeom.png" />

The geometry for the tutorial is taken as a simple wedge with three different boundaries:
* wedge,
* bottom, and
* side.

### Boundary conditions
For the wedge boundary, we prescribe standard OpenFOAM _wedge_ boundary condition for all the variables. For the mass transfer (`omegaV` and `omegaC` variables), we prepared custom Robin external mass transfer boundary conditions _breadOmegaVMixed_. The boundary conditions can be changed similarly as in other OpenFOAM software in `0.org/` directory. The boundary condition for the bottom, side and top patches differ only in the external mass transfer coefficient `kM`. Furthermore, the temporal evolution of the water vapors and carbon dioxide in oven can be changed in `constant/omegaVInfTable` `constant/omegaCInfTable` as standard OpenFOAM interpolation table.

`0.org/omegaV`
```
boundaryField
{
    // Zhang experiment
    "(wedgeZ0|wedgeZE)"
	{
		type wedge;
	}   
    "symmetryPatch"
    {
        type  symmetry;
    }
    
    sides
    {
        type breadOmegaVMixed;;
        kM               1e-5; // -- mass transfer coefficient

        // -- mixed BC mandatory entires
        refValue        uniform 1e-3;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 1e-3;
        omegaVInfTableDict   
        {
            file "$FOAM_CASE/constant/omegaCInfTable";
            outOfBounds warn;
        }
        DFieldName      DEffcM; // -- name of the diffusion coeficient field
	}

    bottom
    {
        type breadOmegaVMixed;;
        kM               0.01; // -- mass transfer coefficient

        // -- mixed BC mandatory entires
        refValue        uniform 1e-3;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 1e-3;
        omegaVInfTableDict   
        {
            file "$FOAM_CASE/constant/omegaCInfTable";
            outOfBounds warn;
        }
        DFieldName      DEffcM; // -- name of the diffusion coeficient field
	}
}
```

`0.org/omegaC`
```
boundaryField
{
    // Zhang experiment
    "(wedgeZ0|wedgeZE)"
	{
		type wedge;
	}   
    "symmetryPatch"
    {
        type  symmetry;
    }
    
    sides
    {
        type breadOmegaVMixed;;
        kM               0.01; // -- mass transfer coefficient

        // -- mixed BC mandatory entires
        refValue        uniform 1e-3;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 1e-3;
        omegaVInfTableDict   
        {
            file "$FOAM_CASE/constant/omegaCInfTable";
            outOfBounds warn;
        }
        DFieldName      DEffcM; // -- name of the diffusion coeficient field
	}

    bottom
    {
        type breadOmegaVMixed;;
        kM               0.01; // -- mass transfer coefficient

        // -- mixed BC mandatory entires
        refValue        uniform 1e-3;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 1e-3;
        omegaVInfTableDict   
        {
            file "$FOAM_CASE/constant/omegaCInfTable";
            outOfBounds warn;
        }
        DFieldName      DEffcM; // -- name of the diffusion coeficient field
	}
}
```

for pressure fixed Dirichlet value is prescribed

`0.org/pG`

```
boundaryField
{
    "symmetryPatch"
    {
        type  symmetry;
    }
    
    sides
    {
        type fixedValue;
        value $internalField;
    }
    "(wedgeZ0|wedgeZE)"
	{
		type wedge;
	}
    bottom
    {
        type fixedValue;
        value $internalField;
    }
}
```
Similarly for the temperature, the external transport in the oven is approximated by the custom Robin boundary condition which is assumed to be same at all boundaries and can be changed in `0.org/T`.
``` 
boundaryField
{
	
    "symmetryPatch"
    {
        type  symmetry;
    }
    sides
	{
        type breadTMixed;
        refValue        uniform 300;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 300;
        alpha           10;
        TInfTableDict   
        {
            file "$FOAM_CASE/constant/TInfTable";
            outOfBounds warn;
        }
	}

    bottom
	{
        type breadTMixed;
        // type breadTBottom;
        refValue        uniform 300;
        refGradient     uniform 0;
        valueFraction   uniform 0;
        value           uniform 300;
        alpha           10;
        TInfTableDict   
        {
            file "$FOAM_CASE/constant/TInfTableBottom";
            outOfBounds warn;
        }
	}
    "(wedgeZ0|wedgeZE)"
	{
		type wedge;
	}
}
```

Here, `alpha` is the external heat transfer coefficient, and again, the temporal evolution of the oven temperature (i.e. baking curve) can be set in `constant/TInfTable` as OpenFOAM interpolation table. Finally, the _fixedDisplacementZeroShear_ boundary condition is prescribed for the deformation at the bottom  and side patches, and the custom _breadDSide_ boundary condition is prescribed for the top patch.

```
top
{
    type            breadDSide;
    sidePos         5.8e-2;
    refValue        uniform (0 0 0);
    refGradient     uniform (0 0 0);
    valueFraction   uniform 1;
    value           uniform 0;
}
```
This boundary condition acts as solids4foam _solidTraction_ boundary condition with zero traction and pressure, until the `sidePos` in the horizontal dimension is reached by some face. Then, this face does not move anymore.  

Alternatively, you can change external heat and mass transfer coefficients for all the boundaries in `'''Boundary conditions'''` section of the  `pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py` control script.
```
'''Boundary and initial conditions'''
TProofing = 301
TStart = 298
TTop = 200
TBottom = 200

kMSidesOmega = 0.01 # -- external mass transfer coeficient 
kMBottomOmega = 0.01

alphaG = 10 # -- external heat transfer coeficient 
alphaGBottom = 14
```

### Internal transfer parameters
The parameters for the internal transfer in the bread can be changed directly in the `constant/transportProperties` and `constant/thermophysicalProperties` or in `'''Internal transport parameters'''` section of the control python script (`pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py`).

```
    '''Internal transport parameters'''
    DFree = 2.6e-5    # -- free volumetric difusivity of the water vapors in CO2 at 300 K
    Dl = 6e-11  # -- liquid water difusivity in the dough
    tortOpen = 2.4   # -- tortuosity
    tortClosed = 70   # -- tortuosity (not used)

    # -- heat conductivity of the dough material with porosity 0, i.e. the 
    # -- absolute term in equation (5) in 
    # -- https://doi.org/10.1016/j.fbp.2008.04.002
    lambdaS = 0.42  # -- heat conductivity of the solid phase (works with addiditional)

    # -- closed-cell bread intristic permeability
    perm = 0.9e-15  # -- bread permeability 

    # -- heat capacities for the individual phases
    CpS = 1130   # -- solid phase
    CpG = 853  # -- CO2
    CpVapor = 1878 # -- water vapors
    CpL = 4200  # -- liquid phase

    # -- mass density for the individual phases
    rhoS = 865

    # -- initial dough volumetric fraction
    alphaD0 = 0.91 
    if not proofing:
        alphaD0 = 0.43
        rhoS = 1387
```

`DFree` parameter sets up the free volumetric diffusivity of the water vapors in carbon dioxide. The temperature and composition dependence of the effective diffusivity is then calculated directly in the solver. `lambdaS` sets up the heat conductivity of the dough material with zero porosity, i.e. the absolute term in equation (5) in https://doi.org/10.1016/j.fbp.2008.04.002 that is used for calculation of the effective heat conductivity. Specific heat capacities and mass densities can be then changed by `Cp` and `rho` parameters.

### Evaporation and fermentation
Evaporation is calculated using Hertz-Knudsen equation while the needed water activity is calculated using Oswin model with parameters measured in https://doi.org/10.1016/0260-8774(91)90020-S. Fermentation kinetics is taken directly from equation (32) in https://doi.org/10.1002/aic.10518. The parameters for all the relations for evaporation and fermentation evaluation can be changed in `constant/reactiveProperties` file or in `'''Evaporation and CO2 generation parameters'''` section of the control python script (`pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py`).
```
    '''Evaporation and CO2 generation parameters'''
    # -- evaporation / condensation coeficient in Hertz-Knudsen equation
    kMPCOpen = 0.015
    kMPCClosed = 0.015

    # -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S) (legacy -- not used)
    evCoef1 = -0.0071
    evCoef2 = 4.5
    n = 0.38

    # -- pre-exponential factor and Tm in CO2 generation kinetics in equation (32) in https://doi.org/10.1002/aic.10518  in (kg/m3/s)
    R0 = 2.3e-3  
    Tm = 313
    deltaT = 14
```
`kMPC` sets up the evaporation coefficient in the Hertz-Knudsen formula. `evCoef1` and `evCoef` are the coefficients for the Oswin model for water activity. Finally, `R0` and `Tm` are the pre-exponential factor and temperature of the fermentation maximum in CO2 generation kinetics.

### Mechanical properties
Bread equalibrium shear and viscous moduli, Poisson ratio and relaxation time can be changed directly in `constant/mechanicalProperties` file or in `'''Mechanical properties'''` section of the control python script (`pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py`).
```
    '''Mechanical properties'''
    withDeformation = 1 # -- turn on (1) /off (0) deformation

    nu = 0.14   # -- Poisson ratio
    E = 30000   # -- Youngs modulus (legacy -- not used) 
    mu0Raw = 147   # kappa = 2*mu*nu/(1-2*nu)   
    muV1Raw = 7400 
    bakedCoeff = 25
    tau1 = 1.6
```

## Running the tutorial
As written above, the tutorial can be either run directly by `Allrun` script in tutorial directory `tutorials/breadAx2DOurExp` or by control python script `pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py` which allows further setup. 
```
# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/breadAx2DOurExp/' # -- base case for simulation
outFolder = '../ZZ_cases/01_breadAx2DOurExp/V9_DlTDep/'

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
runPostProcess = True   # -- run post-processing

proofing = False  # -- proofing included
proofing = True  # -- proofing included
```
`baseCaseDir` sets up the tutorial directory, `outFolder` specifies path where the tutorial will be copied, modified and run. 

### Parallel run
The tutorial is prepared to run also in parallel. It is possible to run it by changing `nCores` parameter in `pyCtrlScripts/runBread2DOurExpFreeBreadProofing.py` to number higher than 1.

