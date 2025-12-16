#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python script to set up and run bread baking simulations according to 
# Zhang et al. https://doi.org/10.1002/aic.10518

# IMPORTS===============================================================
import sys
from OF_caseClass import OpenFOAMCase
import numpy as np
from blockMeshDictClassV8 import *
from meshGeneration import *
import re
import matplotlib.pyplot as plt

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/breadAx2D/' # -- base case for simulation
outFolder = '../ZZ_cases/00_breads/breadAx2D/'

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
runPostProcess = True   # -- run post-processing

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.1e-2 # -- aproximate computational cell size
rLoaf = 3.6e-2  # -- loaf radius                
hLoaf = 3.5e-2  # -- loaf height
arcL = 0.008    # -- length of the arc at the side of the bread   

'''Internal transport parameters'''
# -- free volumetric difusivity of the water vapors in CO2 at 300 K
DFree = 2.2e-6 

# -- heat conductivity of the dough material with porosity 0, i.e. the 
# -- absolute term in equation (5) in 
# -- https://doi.org/10.1016/j.fbp.2008.04.002
lambdaS = 0.45 

perm = 0.9e-12  # -- bread permeability 

# -- heat capacities for the individual phases
CpS = 1200   # -- solid phase
CpG = 853  # -- CO2
CpVapor = 1878 # -- water vapors
CpL = 4200  # -- liquid phase

# -- mass densities for the individual phases
rhoS = 700  # -- solid density    
rhoL = 1000  # -- liquid density   

'''Evaporation and CO2 generation parameters'''
# -- evaporation / condensation coeficient in Hertz-Knudsen equation
kMPC = 0.03

# -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S)
evCoef1 = -0.0056
evCoef2 = 5.5

# -- pre-exponential factor and Tm in CO2 generation kinetics 
# -- in equation (32) in https://doi.org/10.1002/aic.10518
R0 = 23e-4 
Tm = 314

'''Mechanical properties'''
withDeformation = 1 # -- turn on (1) /off (0) deformation
nu = 0.15   # -- Poisson ratio
E = 12000   # -- Youngs modulus

'''Numerics'''
timeStep = 1    # -- computational time step
plusTime1 = 880 # -- how long to run with deformation
plusTime2 = 20 # -- how long to run without deformation
writeInt = 20   # -- how often to write results
nIter = 50  # -- number of iterations in each time step
dynSolver = 'breadBakingFoam'   # -- used solver
nCores = 4 # -- number of cores to run the simulation

# -- relaxation factors
DRelax = 0.1
DFinalRelax = 1

'''Boundary conditions'''
kG = 0.01   # -- external mass transfer coeficient
alphaG = 10 # -- external heat transfer coeficient 

'''Post-processing'''
fig, axs = plt.subplots(3, 1, figsize=(9, 16))  # figure with plots

# SCRIPT ITSELF (DO NOT EDIT)===========================================                       
# -- create OpenFOAMCase object to change values in dictionaries
baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase(baseCaseDir)
baseCase.changeOFCaseDir(outFolder)
baseCase.copyBaseCase()

# OTHER COMPUTATIONS====================================================
dA = mSStep
dX, dY, dZ = dA, dA, dA                                  
x0 = y0 = z0 = 0.0      
grX = grY = grZ = "1.0"

# -- prepare blockMeshDict using luckas python class
if prepBlockMesh:
    prep2DMeshZhang(arcL, rLoaf, hLoaf, x0, y0, z0, dA, dX, dY, dZ, grX, grY, grZ, baseCase)

# CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
# 1) BOUNDARY CONDITIONS
# -- change in tutorial case

# 2) constant/transportProperties
baseCase.setParameters(
    [
        ['constant/transportProperties', 'withDeformation', str(withDeformation), ''],
        ['constant/transportProperties', 'permGLViscG', str(perm), ''],
    ]
)

# 3) constant/thermophysicalProperties
baseCase.setParameters(
    [
        ['constant/thermophysicalProperties', 'lambda', str(lambdaS), 'solid'],
        ['constant/thermophysicalProperties', 'rho', str(rhoS), 'solid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpS), 'solid'],
        ['constant/thermophysicalProperties', 'rho', str(rhoL), 'liquid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpL), 'liquid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpG), 'CO2'],
        ['constant/thermophysicalProperties', 'Cp', str(CpVapor), 'vapor'],
        ['constant/thermophysicalProperties', 'D', str(DFree), 'transport'],
    ]
)

# 4) constant/reactiveProperties
# -- parameters for evaporation and CO2 generation
baseCase.setParameters(
    [
        ['constant/reactiveProperties', 'kMPCOpen', str(kMPC), 'evaporation'],
        ['constant/reactiveProperties', 'kMPCClosed', str(kMPC), 'evaporation'],
        ['constant/reactiveProperties', 'evCoef1', str(evCoef1), 'evaporation'],
        ['constant/reactiveProperties', 'evCoef2', str(evCoef2), 'evaporation'],
        ['constant/reactiveProperties', 'R0', str(R0), 'fermentation'],
        ['constant/reactiveProperties', 'Tm', str(Tm), 'fermentation'],
    ]
)
        
# 5 system/controlDict
baseCase.setParameters(
    [
        ['system/controlDict', 'endTime', str(plusTime1), ''],
        ['system/controlDict', 'deltaT', '%.5g'%timeStep, ''],
        ['system/controlDict', 'writeInterval', '%.5g'%writeInt, ''],
    ]
)

# 6) fvSolutions
baseCase.setParameters(
    [
        ['system/fvSolution', 'nOuterCorrectors', str(nIter), 'PIMPLE'],
        ['system/fvSolution', 'D', str(DRelax), 'fields'],
        ['system/fvSolution', 'DFinal', str(DFinalRelax), 'fields'],
    ]
)

# 7) mechanical properties
baseCase.setParameters(
    [
        ['constant/mechanicalProperties', 'nu', str(nu), 'bread'],
        ['constant/mechanicalProperties', 'E', str(E), 'bread']
    ]
)

# -- prepare geom
if makeGeom:
    baseCase.runCommands(
        [
            'chmod 755 ./* -R',
            'blockMesh > log.blockMesh',
            'rm -rf 0',
            'cp -r 0.org 0',
            'paraFoam -touch',
        ]
    )

# RUN THE SIMULATION====================================================
if runDynSim:
    if nCores > 1:
        baseCase.setParameters(
            [
                ['system/decomposeParDict', 'numberOfSubdomains', str(nCores), '']
            ]
        )
        baseCase.runCommands(
            [
                'decomposePar > log.decomposePar',
                'foamJob -parallel -screen %s > log.%s' %(dynSolver,dynSolver),
            ]
        )
    else:
        baseCase.runCommands(
            [
                '%s > log.%s' %(dynSolver,dynSolver),
            ]
        )

    # -- run the rest of the simualation without further deformation
    if plusTime2 > 0:
        baseCase.setParameters(
            [
                ['system/controlDict', 'endTime', str(plusTime1 + plusTime2), ''],
                ['constant/transportProperties', 'withDeformation', '0', '']
            ]
        )
        if nCores > 1:
            baseCase.runCommands(
                [
                    'foamJob -parallel -screen %s > log.%s_2' %(dynSolver,dynSolver),
                ]
            )
        else:
            baseCase.runCommands(
                [
                    '%s > log.%s_2' %(dynSolver,dynSolver),
                ]
            )
        
# POST-PROCESSING=======================================================
if runPostProcess:
    # -- load the experimental data
    TExpCenter = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_Zhang_center.dat', delimiter=';')
    TExpSurface = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_Zhang_surface.dat', delimiter=';')
    DExp = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_Zhang_D.dat', skiprows=1)
    moistureExp = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_Zhang_moisture.dat', skiprows=1, delimiter=';')
    
    # -- run post-processing tasks
    if nCores == 1:
        baseCase.updateTimes()
        baseCase.runCommands(
            [
                'postProcess -func "probeZhang" -dict system/probeZhang > log.postProcess',
                'rm -rf 0',
                'intMoisture > log.intMoisture',
            ]
        )
    else:
        baseCase.updateTimesParallel()
        baseCase.runCommands(
            [
                'foamJob -parallel -screen postProcess -func "probeZhang" -dict system/probeZhang > log.postProcess',
                'rm -rf processor*/0',
                'foamJob -parallel -screen intMoisture > log.intMoisture',
            ]
        )

    # -- gather the displacement data from probe points
    rows = []
    lines = []
    D = []
    nProbes = 3
    if nCores > 1:
        latestTime = baseCase.latestParTime
    else:
        latestTime  = baseCase.latestTime
    with open(baseCase.dir + '/postProcessing/probeZhang/%d/D'%latestTime, 'r') as fl:
        lines = fl.readlines()
        lines = lines[nProbes+2:]
        # print(lines)
        for line in lines:
            parts = line.split(") (")
            first_entry = parts[0].split(maxsplit=1)
            vectors = [first_entry[1]] if len(first_entry) > 1 else []
            vectors.extend(parts[1:])

            vectors = [
                tuple(map(float, vec.replace("(", "").replace(")", "").split()))
                for vec in vectors
            ]
            rows.append(vectors)

    # -- Convert displacements to numpy array
    D = np.array(rows)
        
    # -- Load temperature profiles in probe points
    probesT = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/T'%latestTime, skiprows=3)

    # -- Load total moisture evolution 
    file_path = "%s/log.intMoisture" %baseCase.dir
    skiprows = -1
    endLine = -1
    with open(file_path, "r") as file:
        lines = file.readlines()
    for lineI in range(len(lines)):
        if 'Time = ' in lines[lineI]:
            skiprows = lineI
            break
    for lineI in range(len(lines)):
        if 'End' in lines[lineI]:
            endLine = lineI
            break

    # Extract all numbers using regex
    numbers = []
    for lineI in range(skiprows, endLine):
        line = lines[lineI]
        matches = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)  # Matches integers and decimals
        numbers.extend(map(float, matches))  # Convert to float and add to the list

    # Convert the list to a NumPy array
    moistureSim = np.array(numbers).reshape(-1,2)

    # -- Temperatures
    axs[0].plot(TExpCenter[:,0],TExpCenter[:,1], 'xr',  label='center temperature experiment')
    axs[0].plot(TExpSurface[:,0],TExpSurface[:,1], 'xb', label='surface temperature experiment')
    axs[0].plot(probesT[:,0] / 60, probesT[:,1] - 273, 'r', label='center temperature simulation')
    axs[0].plot(probesT[:,0] / 60, probesT[:,2] - 273, 'b', label='center temperature simulation')
    axs[0].set_xlabel("time (min)")
    axs[0].set_ylabel("T (°C)")
    axs[0].set_xlim(0, 15)
    axs[0].set_title("Temperature evolution in the center and at the surface")
    axs[0].legend()

    # -- Moisture
    axs[1].plot(moistureSim[:,0] / 60, moistureSim[:,1], 'b', label='simulation')
    axs[1].plot(moistureExp[:,0] , moistureExp[:,1], 'xb', label='experiment')
    axs[1].set_xlabel("time (min)")
    axs[1].set_ylabel("total moisture content (-)")
    axs[1].set_xlim(0,15)
    axs[1].set_title("Total moisture content in the the bread")
    axs[1].legend()

    # -- Displacement
    axs[2].plot(probesT[1:,0] / 60, D[:, 1, 1], 'b', label='simulation DY')
    axs[2].plot(probesT[1:,0] / 60, D[:, 2, 0], 'r', label='simulation DX')
    axs[2].plot(DExp[:,0] / 60, DExp[:,1], 'xr', label='experimental DX')
    axs[2].plot(DExp[:,0] / 60, DExp[:,2], 'xb', label='experimental DY')
    axs[2].set_xlabel("time (min)")
    axs[2].set_ylabel("displacement in X and Y directions")
    axs[2].set_xlim(0,15)
    axs[2].set_title("Displecement of the bread in vertical (X) and horizontal (Y) directions")
    axs[2].legend()
    fig.tight_layout()

    plt.savefig(baseCase.dir + 'postProcessingPlot.png')
                                        
