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
from myAddFcs import *
import re
import matplotlib.pyplot as plt

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/bread3DOurExp/' # -- base case for simulation
# baseCaseDir = '../ZZ_cases/00_breads/breadAx2DOurExp/'
# baseCaseDir = '../ZZ_cases/00_breads/fine_breadAx2DOurExp_lam0.44/'
# baseCaseDir = '../ZZ_cases/00_breads/bread3DOurExp/'
baseCaseDir = '../ZZ_cases/00_breads/bread2DOurFree/'
outFolder = '../ZZ_cases/00_breads/bread2DOurFree/'

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
prepBlockMesh = False    # -- preparation of the blockMeshDict script
makeGeom = False # -- creation of the geometry for computation
runDynSim = False    # -- run simulation
runPostProcess = True   # -- run post-processing
# runPostProcess = False   # -- run post-processing

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.1e-2 # -- aproximate computational cell size
rLoaf1 = 8.0e-2  # -- loaf radius                
rLoaf2 = 8.0e-2  # -- loaf radius                
hLoaf = 7e-2  # -- loaf height 

'''Internal transport parameters'''
# -- free volumetric difusivity of the water vapors in CO2 at 300 K
DFree = 2.22e-5 

# -- heat conductivity of the dough material with porosity 0, i.e. the 
# -- absolute term in equation (5) in 
# -- https://doi.org/10.1016/j.fbp.2008.04.002
lambdaS = 0.44

perm = 3e-12  # -- bread permeability 

# -- heat capacities for the individual phases
CpS = 700   # -- solid phase
CpG = 853  # -- CO2
CpVapor = 1878 # -- water vapors
CpL = 4200  # -- liquid phase

# -- mass densities for the individual phases
rhoS = 700  # -- solid density    
rhoL = 1000  # -- liquid density   

'''Evaporation and CO2 generation parameters'''
# -- evaporation / condensation coeficient in Hertz-Knudsen equation
kMPC = 0.023

# -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S)
evCoef1 = -0.0056
evCoef2 = 5.5

# -- pre-exponential factor and Tm in CO2 generation kinetics 
# -- in equation (32) in https://doi.org/10.1002/aic.10518
R0 = 6e-4 
Tm = 314
Tm = 308

'''Mechanical properties'''
withDeformation = 1 # -- turn on (1) /off (0) deformation
nu = 0.15   # -- Poisson ratio
E = 12000   # -- Youngs modulus

'''Numerics'''
timeStep = 0.5  # -- computational time step
plusTime1 = 870 # -- how long to run with deformation
plusTime2 = 730 # -- how long to run without deformation
writeInt = 10   # -- how often to write results
nIter = 50  # -- number of iterations in each time step
dynSolver = 'breadBakingFoam'   # -- used solver
nCores = 4 # -- number of cores to run the simulation

# -- relaxation factors
DRelax = 0.2
DFinalRelax = 1

'''Boundary conditions'''
kMSides = 6e-4   # -- external mass transfer coeficient
kMBottom = 3e-4   # -- external mass transfer coeficient
kMTop = 0.01   # -- external mass transfer coeficient
alphaG = 10 # -- external heat transfer coeficient 

'''Post-processing'''
fig, axs = plt.subplots(3, 2, figsize=(18, 16))  # figure with plots

# SCRIPT ITSELF (DO NOT EDIT)===========================================                       
# -- create OpenFOAMCase object to change values in dictionaries
baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase(baseCaseDir)
baseCase.changeOFCaseDir(outFolder)
# baseCase.copyBaseCase()

# OTHER COMPUTATIONS====================================================
dA = mSStep
dX, dY, dZ = dA, dA, dA                                  
x0 = y0 = z0 = 0.0      
grX = grY = grZ = "1.0"

# -- prepare blockMeshDict using luckas python class
if prepBlockMesh:
    prep3DMeshOurExp(rLoaf1, rLoaf2, hLoaf, dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True)

# CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
# 1) BOUNDARY CONDITIONS
# -- change in tutorial case
# baseCase.setParameters(
#     [
#         ['0.org/omegaV', 'kM', str(kMSides), 'sides'],
#         ['0.org/omegaV', 'kM', str(kMBottom), 'bottom'],
#         ['0.org/omegaV', 'kM', str(kMTop), 'top'],
#         ['0.org/pG', 'kM', str(kMSides), 'sides'],
#         ['0.org/pG', 'kM', str(kMBottom), 'bottom'],
#     ]
# )

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
            'extrudeMesh > log.extrudeMesh',
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
    expData = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_all.dat', skiprows=1)
    expData2 = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_all_2.dat', skiprows=1)
    
    # -- run post-processing tasks
    if nCores == 1:
        baseCase.updateTimes()
        baseCase.runCommands(
            [
                'postProcess -func "probeOur" -dict system/probeOur > log.postProcess',
                'TLFProbe -point "(0.013 0.041 0)" > log.TPoint1',
                'TLFProbe -point "(0.035 1e-4 0)" > log.TPoint2',
                'TLFProbe -point "(0.016 1e-4 0)" > log.TPoint3',
                'TLFProbe -point "(0.028 0.021 0)" > log.TPoint4',
                'rm -rf 0',
                'intMoisture > log.intMoisture',
            ]
        )
    else:
        baseCase.updateTimesParallel()
        baseCase.runCommands(
            [
                'rm -rf processor*/0',
                'foamJob -parallel -screen postProcess -func "probeOur" -dict system/probeOur > log.postProcess',
                'foamJob -parallel -screen TLFProbe -point "(0.012 1e-4 0)" > log.TPoint6',
                'foamJob -parallel -screen TLFProbe -point "(0.061 1e-3 0)" > log.TPoint7',
                'foamJob -parallel -screen TLFProbe -point "(0.027 0.047 0)" > log.TPoint5',
                'foamJob -parallel -screen TLFProbe -point "(0.032 0.041 0)" > log.TPoint8',
                'foamJob -parallel -screen intMoisture > log.intMoisture',
            ]
        )

    # -- gather the displacement data from probe points
    rows = []
    lines = []
    D = []
    nProbes = 1
    if nCores > 1:
        latestTime = baseCase.latestParTime
    else:
        latestTime  = baseCase.latestTime
    with open(baseCase.dir + '/postProcessing/probeOur/%d/D'%latestTime, 'r') as fl:
        lines = fl.readlines()
        lines = lines[nProbes+1:]
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
    # probesT = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/T'%latestTime, skiprows=3)

    # -- Load total moisture evolution 
    TPoint6 = readDataFromLogFile("%s/log.TPoint6" %baseCase.dir)
    TPoint7 = readDataFromLogFile("%s/log.TPoint7" %baseCase.dir)
    TPoint5 = readDataFromLogFile("%s/log.TPoint5" %baseCase.dir)
    TPoint8 = readDataFromLogFile("%s/log.TPoint8" %baseCase.dir)
    moistureSim = readDataFromLogFile("%s/log.intMoisture" %baseCase.dir)


    # -- Temperatures
    axs[0,0].plot(expData[:,-2],expData[:,0], '--r',  label='exp. point 5')
    axs[0,0].plot(expData[:,-2],expData[:,1], '--g',  label='exp. point 6')
    axs[0,0].plot(expData[:,-2],expData[:,2], '--b',  label='exp. point 7')
    axs[0,0].plot(expData[:,-2],expData[:,3], '--m',  label='exp. point 8')
    axs[0,1].plot(expData2[:,-2],expData2[:,0], '--r',  label='exp. point 5')
    axs[0,1].plot(expData2[:,-2],expData2[:,1], '--g',  label='exp. point 6')
    # axs[0,1].plot(expData2[:,-2],expData2[:,2], '--b',  label='exp. point 7')
    # axs[0,1].plot(expData2[:,-2],expData2[:,3], '--m',  label='exp. point 8')
    # axs[0].plot(TExpPoint2[:,-1],TExpPoint2[:,0], '--b',  label='exp. point 2')
    # axs[0].plot(TExpPoint3[:,-1],TExpPoint3[:,0], '--g',  label='exp. point 4')
    # axs[0].plot(TExpSurface[:,0],TExpSurface[:,1], 'xb', label='surface temperature experiment')
    axs[0,0].plot(TPoint5[:,0] / 60, TPoint5[:,1] - 273, 'r', label='sim. point 5')
    axs[0,0].plot(TPoint6[:,0] / 60, TPoint6[:,1] - 273, 'g', label='sim. point 6')
    axs[0,0].plot(TPoint7[:,0] / 60, TPoint7[:,1] - 273, 'b', label='sim. point 7')
    axs[0,0].plot(TPoint8[:,0] / 60, TPoint8[:,1] - 273, 'm', label='sim. point 8')
    axs[0,1].plot(TPoint5[:,0] / 60, TPoint5[:,1] - 273, 'r', label='sim. point 5')
    axs[0,1].plot(TPoint6[:,0] / 60, TPoint6[:,1] - 273, 'g', label='sim. point 6')
    axs[0,1].plot(TPoint7[:,0] / 60, TPoint7[:,1] - 273, 'b', label='sim. point 7')
    axs[0,1].plot(TPoint8[:,0] / 60, TPoint8[:,1] - 273, 'm', label='sim. point 8')
    # axs[0].plot(probesT[:,0] / 60, probesT[:,2] - 273, 'b', label='center temperature simulation')
    axs[0,0].set_xlabel("time (min)")
    axs[0,0].set_ylabel("T (°C)")
    axs[0,0].set_ylim(20,120)
    axs[0,0].set_xlim(0, 35)
    axs[0,0].set_title("Temperature evolution in the center and at the surface")

    axs[0,1].set_xlabel("time (min)")
    axs[0,1].set_ylabel("T (°C)")
    axs[0,1].set_ylim(20,120)
    axs[0,1].set_xlim(0, 35)
    axs[0,1].set_title("Temperature evolution in the center and at the surface")
    axs[0,0].legend()

    # -- Moisture
    axs[1,0].plot(moistureSim[:,0] / 60, moistureSim[:,1], 'b', label='simulation')
    axs[1,0].plot(expData[:,-2], expData[:,-1], '--b', label='experiment')
    axs[1,0].set_xlabel("time (min)")
    axs[1,0].set_ylabel("total moisture content (-)")
    axs[1,0].set_ylim(0.35,0.56)
    # axs[1].set_xlim(0,28)
    axs[1,0].set_title("Total moisture content in the the bread")
    axs[1,0].legend()

    axs[1,1].plot(moistureSim[:,0] / 60, moistureSim[:,1], 'b', label='simulation')
    axs[1,1].plot(expData2[:,-2], expData2[:,-1], '--b', label='experiment')
    axs[1,1].set_xlabel("time (min)")
    axs[1,1].set_ylabel("total moisture content (-)")
    axs[1,1].set_ylim(0.35,0.56)
    # axs[1].set_xlim(0,28)
    axs[1,1].set_title("Total moisture content in the the bread")
    axs[1,1].legend()

    # -- Displacement
    axs[2,0].plot(TPoint5[:,0] / 60, D[:, 0, 0], 'b', label='simulation DX')
    # axs[2].plot(DExp[:,0] / 60, DExp[:,1]+1.75e-2, 'xb', label='experimental DX')
    # axs[2].plot(DExp[:,0] / 60, DExp[:,2], 'xb', label='experimental DY')
    axs[2,0].set_xlabel("time (min)")
    axs[2,0].set_ylabel("displacement in X and Y directions")
    axs[2,0].set_xlim(0,35)
    axs[2,0].set_title("Displecement of the bread in vertical (X) and horizontal (Y) directions")
    axs[2,0].legend()
    fig.tight_layout()

    plt.savefig(baseCase.dir + 'postProcessingPlot.png')
                                        
