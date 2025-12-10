# 2D axisymmetrical bread acording to Zhang
## Case description and setup
This tutorial shows a two-dimensional internal simulation of the bread in the oven. External transport is resolved by custom mixed boundary conditions. The tutorial is located in `tutorials/breadAx2D` and can be:
1. run directly as prepared by `Allrun` script in `tutorials/breadAx2D` folder, or
2. modified and run by `pyCtrlScripts/runBreadAx2D.py` control script.

The description of the solved equations and variables is in greater detail discussed in https://doi.org/10.14311/TPFM.2025.015. Furthermore in the solver, the solved variables are noted as: 
* `alphaI` - volumetric fraction of the I-th phase,
* `T` - temperature,
* `pG` - pressure of the gas phase,
* `omegaV` - mass fraction of the water vapors in the gas, and
* `D` - deformation vector. 

### Geometry and computational mesh description

<img alt="twoBreadsLength" src="twoBreadsLengthV2.png" />

Geometry for the tutorial is taken from the work of Zhang (https://doi.org/10.1002/aic.10518). The computational mesh is prepared using `system/blockMeshDict`, which is pre-prepared in the tutorial case. `blockMeshDict` file with the modified dimesions (radius, height and lenght of the arc) or the computational cell size can be generated using `pyCtrlScripts/runBreadAx2D.py` script by changing `'''Geometry parameters'''` part:
```
'''Geometry parameters'''
typeOfMesh = '2DZhang'
mSStep = 0.1e-2 # -- aproximate computational cell size
rLoaf = 3.6e-2  # -- loaf radius                
hLoaf = 3.5e-2  # -- loaf height
arcL = 0.008    # -- length of the arc at the side of the bread
```

There are three diferent types of the geometry boundaries:
* wedge (depicted in green),
* bottom (depicted in blue), and 
* side (depicted in red).

### Boundary conditions
For the wedge boundary, we prescribe standard OpenFOAM _wedge_ boundary condition for all the variables. The bread in the experiment is placed on the metal grid which allows the external mass transfer also from the bottom side. Therefore, all the variables but the deformation are prescribed with the same boundary condition for both the side and bottom patches in this tutorial. For the mass transfer (`pG` and `omegaV` variables), we prepared custom Robin external mass transfer boundary conditions _breadPGMixed_ and _breadOmegaVMixed_, respectively. The boundary conditions can be changed similarly as in other OpenFOAM software in `0.org/` directory. Note that for _pG_, the Dirichlet boundary condition `pG` = 100 000 Pa at both _bottom_ and _side_ patches is used. For the mass fraction of the water vapors in the gas _omegaV_, the external mass transfer coefficient `kM` can be changed in `0.org/omegaV`

```
"(sides|bottom)"
{
    type breadOmegaVMixed;
    kM               0.01; // -- mass transfer coefficient

    // -- mixed BC mandatory entires
    refValue        uniform 7.7e-3;
    refGradient     uniform 0;
    valueFraction   uniform 0;
    value           uniform 7.7e-3;
    omegaVInfTableDict   
    {
        file "$FOAM_CASE/constant/omegaVInfTable";
        outOfBounds warn;
    }
}
```
Next, the temporal evolution of the water vapors in the oven can be changed in `constant/omegaVInfTable` using standard OpenFOAM interpolation table. Similarly for the temperature, the external transport in the oven is approximated by the custom Robin boundary condition which can be changed in `0.org/T`
``` 
"(sides|bottom)"
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
```
Here, `alpha` is the external heat transfer coefficient, and again, the temporal evolution of the oven temperature (i.e. baking curve) can be set in `constant/TInfTable` as OpenFOAM interpolation table. Finally, the Dirichlet zero boundary condition is prescribed for the deformation and the bottom patch, and the custom _breadDFloor_ boundary condition is prescribed for the side patch.
```
"(sides)"
{
    type            breadDFloor;
    floorPos        -17.5e-3;
    refValue        uniform (0 0 0);
    refGradient     uniform (0 0 0);
    valueFraction   uniform 1;
    value           uniform 0;
}
```
This boundary condition acts as solids4foam _solidTraction_ boundary condition with zero traction and pressure, until the `floorPos` in the vertical dimension is reached by some face. Then, this face does not move anymore.      

## Internal transfer parameters
The parameters for the internal transfer in the bread can be changed directly in the `constant/transportProperties` and `constant/thermophysicalProperties` or in `'''Internal transport parameters'''` section of the control python script (`pyCtrlScripts/runBreadAx2D.py`).

```
'''Internal transport parameters'''
# -- free volumetric difusivity of the water vapors in CO2 at 300 K
DFree = 2.22e-6 

# -- heat conductivity of the dough material with porosity 0, i.e. the 
# -- absolute term in equation (5) in 
# -- https://doi.org/10.1016/j.fbp.2008.04.002
lambdaS = 0.447 

perm = 0.9e-12  # -- bread permeability 

# -- heat capacities for the individual phases
CpS = 700   # -- solid phase
CpG = 853  # -- CO2
CpVapor = 1878 # -- water vapors
CpL = 4200  # -- liquid phase

# -- mass densities for the individual phases
rhoS = 700  # -- solid density    
rhoL = 1000  # -- liquid density   
```

`DFree` parameter sets up the free volumetric diffusivity of the water vapors in carbon dioxide. The temperature and composition dependence of the effective diffusivity is then calculated directly in the solver. `lambdaS` sets up the heat conductivity of the dough material with zero porosity, i.e. the absolute term in equation (5) in https://doi.org/10.1016/j.fbp.2008.04.002 that is used for calculation of the effective heat conductivity. Specific heat capacities and mass densities can be then changed by `Cp` and `rho` parameters.

### Evaporation and fermentation
Evaporation is calculated using Hertz-Knudsen equation while the needed water activity is calculated using Oswin model with parameters measured in https://doi.org/10.1016/0260-8774(91)90020-S. Fermentation kinetics is taken directly from equation (32) in https://doi.org/10.1002/aic.10518. The parameters for all the relations for evaporation and fermentation evaluation can be changed in `constant/reactiveProperties` file or in `'''Evaporation and CO2 generation parameters'''` section of the control python script (`pyCtrlScripts/runBreadAx2D.py`).
```
'''Evaporation and CO2 generation parameters'''
# -- evaporation / condensation coeficient in Hertz-Knudsen equation
kMPC = 0.42

# -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S)
evCoef1 = -0.0056
evCoef2 = 5.5

# -- pre-exponential factor and Tm in CO2 generation kinetics 
# -- in equation (32) in https://doi.org/10.1002/aic.10518
R0 = 22e-4 
Tm = 314
```
`kMPC` sets up the evaporation coefficient in the Hertz-Knudsen formula. `evCoef1` and `evCoef` are the coefficients for the Oswin model for water activity. Finally, `R0` and `Tm` are the pre-exponential factor and temperature of the fermentation maximum in CO2 generation kinetics.

### Mechanical properties
Bread Youngs modulus and Poisson ratio can be changed directly in `constant/mechanicalProperties` file or in `'''Mechanical properties'''` section of the control python script (`pyCtrlScripts/runBreadAx2D.py`).
```
'''Mechanical properties'''
withDeformation = 1 # -- turn on (1) /off (0) deformation
nu = 0.15   # -- Poisson ratio
E = 12000   # -- Youngs modulus
```

## Running the tutorial
As written above, the tutorial can be either run directly by `Allrun` script in tutorial directory `tutorials/breadAx2D` or by control python script `pyCtrlScripts/runBreadAx2D.py` which allows further setup. 
```
# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/breadAx2D/' # -- base case for simulation
outFolder = '../ZZ_cases/00_breads/breadAx2D/'

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
runPostProcess = True   # -- run post-processing
```
`baseCaseDir` sets up the tutorial directory, `outFolder` specifies path where the tutorial will be copied, modified and run. 

### Parallel run
The tutorial is prepared to run also in parallel. It is possible to run by changing `nCores` parameter in `pyCtrlScripts/runBreadAx2D.py` to number higher than 1.

## Post-processing
### Prepared python post-processing
For the post-processing, it is possible to use prepared `pyCtrlScripts/runBreadAx2D.py` script, which compares the results directly with the experimental data by the generation of the following figure. 

<img alt="tutBreadAx2DPostProcess" src="tutBreadAx2DPostProcess.png" />

The resulting post-processing figure consists of three plots. In the first of them the comparison of the temperature evolution in the center, and at the surface of the bread for the simulation and experiment is depicted. Experimental data are loaded from `ZZ_dataForPostProcessing` directory in the tutorial. In the second plot, the total moisture content in the bread for the simulation and experiment is compared. Finally, in the third plot, the comparison of the simulation and experimental deformation at the bread top (X) and side (Y) are compared.

### Paraview post-processing
To further examine the results, you can visualize them using paraview software. 
1. run the `paraview` (in this description we use Paraview 5.12.1),
2. open `breadAx2D.OpenFOAM` file, that was created during run of `Allrun` script or `pyCtrlScripts/runBreadAx2D.py`,
3. select open data with `Open FOAM Reader`,

<img alt="tutBreadAx2DOpenWith" src="tutBreadAx2DOpenWith.png" />

4. select proper `Case Type` depending on type of your data (single core -> `Reconstructed Case`, parallel -> `Decomposed Case`, 
5. choose the desired `Time` for post-processing (e.g. 240), and
6. select `Apply`.

<img alt="tutBreadAx2DAfterApply" src="tutBreadAx2DAfterApply.png" />

You will obtain such a visualization of the wedge computational mesh. Note that the problem is solved in Total Lagragian formulation, which means mesh is static and the effect of the deformation is resolved by the deformation gradient `F` and its Jacobian `J` in equations.

To vizualize the deformation 
1. press `Ctrl + space` and type `Warp By Vector`, which will apply `Warp By Vector` filter on your data,
2. under `Warp By Vector` filter properties select `pointD` field under `Vectors` to warp the geometry by this field, and click `Apply`,
3. change the view direction to `+Z` and rotate `+90 clockwise`,
4. select 2D interaction mode.

<img alt="tutBreadAx2DAfterWarpByVector" src="tutBreadAx2DAfterWarpByVector.png" />

Now, you can split the layout and vizualize multiple variables at once.

<img alt="tutBreadAx2DAfterVizFields" src="tutBreadAx2DAfterVizFields.png" />

Alternatively, you can try to load prepared paraview state in `ZZ_dataForPostProcessing/seeMultipleFields.pvsm`.
