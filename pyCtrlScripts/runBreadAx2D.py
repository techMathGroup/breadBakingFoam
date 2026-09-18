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
import os

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/breadAx2D/' # -- base case for simulation
# outFolder = '../ZZ_cases/00_breads/newImpl_lambda04const_tortClo40_tortOpen10_3_kmpClosed1_open02_newKin/'

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
# prepBlockMesh = False    # -- preparation of the blockMeshDict script
# makeGeom = False # -- creation of the geometry for computation
# runDynSim = False    # -- run simulation
runPostProcess = True   # -- run post-processing

withKynuti = True
withKynuti = False
timeKynuti = 0
TKynuti = 300
timeStepKynuti = 5

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.1e-2 # -- aproximate computational cell size
rLoaf = 3.6e-2  # -- loaf radius                
hLoaf = 3.5e-2  # -- loaf height
arcL = 0.008    # -- length of the arc at the side of the bread   

'''Internal transport parameters'''
DFree = 2.6e-5    # -- free volumetric difusivity of the water vapors in CO2 at 300 K
Dl = 5e-12  # -- liquid water difusivity in the dough
tortOpen = 2.4   # -- tortuosity
tortClosed = 70   # -- tortuosity

# -- heat conductivity of the dough material with porosity 0, i.e. the 
# -- absolute term in equation (5) in 
# -- https://doi.org/10.1016/j.fbp.2008.04.002
lambdaS = 0.42

# -- intrinsic gas permeabilities of raw and baked bread
# gasPermeabilityRaw = 1.3e-14
# gasPermeabilityRaw = 3.7e-14
gasPermeabilityRaw = 4.5e-14
# gasPermeabilityRaw = 2e-14
gasPermeabilityBaked = 1e-11

# -- heat capacities for the individual phases
CpS = 1130   # -- solid phase
CpG = 853  # -- CO2
CpVapor = 1878 # -- water vapors
CpL = 4200  # -- liquid phase

# -- mass densities for the individual phases
rhoS = 507  # -- solid density    
# rhoS = 764  # -- solid density    
# rhoL = 1000  # -- liquid density   

'''Evaporation and CO2 generation parameters'''
# -- evaporation / condensation coeficient in Hertz-Knudsen equation
kMPCOpen = 0.015
kMPCClosed = 0.015

# -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S) (legacy -- not used)
evCoef1 = -0.0071
evCoef2 = 4.5
n = 0.38

# -- pre-exponential factor and Tm in CO2 generation kinetics in equation (32) in https://doi.org/10.1002/aic.10518 (kg/m^3/s)
R0 = 1.8e-2
# R0 = 3e-2
Tm = 313   
deltaT = 10

# R0 = 1.8e-3
# # Tm = 313
deltaT = 14

'''Mechanical properties'''
withDeformation = 1 # -- turn on (1) /off (0) deformation
# withDeformation = 0 # -- turn on (1) /off (0) deformation
nu = 0.14   # -- Poisson ratio
E = 30000   # -- Youngs modulus
mu0Raw = 230    # kappa = 2*mu*nu/(1-2*nu)  
muV1Raw = 7400
bakedCoeff = 8
tau1 = 1.6


'''Numerics'''
timeStep = 0.25    # -- computational time step
plusTime1 = 360 # -- how long to run with deformation
plusTime2 = 540 # -- how long to run without deformation
writeInt = 10   # -- how often to write results
nIter = 40  # -- number of iterations in each time step
dynSolver = 'breadBakingFoam'   # -- used solver
nCores = 4 # -- number of cores to run the simulation

TRelaxAfter = 0.3

# -- relaxation factors
DRelax = 0.3
DFinalRelax = 1

'''Boundary conditions'''
kG = 0.01   # -- external mass transfer coeficient
alphaG = 10 # -- external heat transfer coeficient 

'''Post-processing'''
fig, axs = plt.subplots(4, 1, figsize=(9, 21))  # figure with plots

outFolder = '../ZZ_cases/00_breads/008_newTau_tau_%g_Dl_%g_kOp_%g_kCl_%g_torOp_%g_torCl_%g_lambdaS_%g_R0_%g_perm_%g/'%(tau1, Dl, kMPCOpen, kMPCClosed, tortOpen, tortClosed, lambdaS, R0, gasPermeabilityRaw)
# baseCaseDir = '../ZZ_cases/00_breads/83_BK15_availSurf_tau_%g_Dl_%g_kOp_%g_kCl_%g_torOp_%g_torCl_%g_lambdaS_%g_R0_%g_perm_%g/'%(tau1, Dl, kMPCOpen, kMPCClosed, tortOpen, tortClosed, lambdaS, R0, gasPermeabilityRaw)


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
kappa0Raw = 2*mu0Raw*nu/(1-2*nu)
muV2Raw = 0
# kappaV1Raw = 2*muV1Raw*nu/(1-2*nu)
kappaV1Raw = 0
mu0Baked = mu0Raw * bakedCoeff
kappa0Baked = 2*mu0Baked*nu/(1-2*nu)
tau2 = 1
tGelat = 65
tau0 = 7

# -- prepare blockMeshDict using luckas python class
if prepBlockMesh:
    prep2DMeshZhang(arcL, rLoaf, hLoaf, x0, y0, z0, dA, dX, dY, dZ, grX, grY, grZ, baseCase)

# CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
# 1) BOUNDARY CONDITIONS
# -- change in tutorial case
baseCase.setParameters(
    [
        ['0.org/T', 'alpha', str(alphaG), 'sides'],
        ['0.org/T', 'alpha', str(alphaG), 'bottom'],
    ]
)

if withKynuti:
    # -- change external temperature
    with open(os.path.join(baseCase.dir, "constant", "TInfTable"), "w") as fl:
        fl.writelines("(\n")
        fl.writelines("\t(0\t%f)\n"%TKynuti)
        fl.writelines("\t(%d\t%f)\n"%(timeKynuti, TKynuti))
        fl.writelines("\t(%f\t%f)\n"%(timeKynuti + 0.1, 463))
        fl.writelines("\t(%d\t%f)\n"%(100000, 463))
        fl.writelines(")\n")

# 2) constant/transportProperties
baseCase.setParameters(
    [
        ['constant/transportProperties', 'withDeformation', str(withDeformation), ''],
        ['constant/transportProperties', 'gasPermeabilityRaw', str(gasPermeabilityRaw), ''],
        ['constant/transportProperties', 'gasPermeabilityBaked', str(gasPermeabilityBaked), ''],
        # ['constant/transportProperties', 'tort', str(tort), ''],
        ['constant/transportProperties', 'tortOpen', str(tortOpen), ''],
        ['constant/transportProperties', 'tortClosed', str(tortClosed), ''],
        
    ]
)

# 3) constant/thermophysicalProperties
baseCase.setParameters(
    [
        ['constant/thermophysicalProperties', 'lambda', str(lambdaS), 'solid'],
        ['constant/thermophysicalProperties', 'rho', str(rhoS), 'solid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpS), 'solid'],
        # ['constant/thermophysicalProperties', 'rho', str(rhoL), 'liquid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpL), 'liquid'],
        ['constant/thermophysicalProperties', 'Cp', str(CpG), 'CO2'],
        ['constant/thermophysicalProperties', 'Cp', str(CpVapor), 'vapor'],
        ['constant/thermophysicalProperties', 'D', str(DFree), 'transport'],
        ['constant/thermophysicalProperties', 'Dl', str(Dl), 'transport'],

    ]
)

# 4) constant/reactiveProperties
# -- parameters for evaporation and CO2 generation
baseCase.setParameters(
    [
        ['constant/reactiveProperties', 'kMPCOpen', str(kMPCOpen), 'evaporation'],
        ['constant/reactiveProperties', 'kMPCClosed', str(kMPCClosed), 'evaporation'],
        ['constant/reactiveProperties', 'evCoef1', str(evCoef1), 'evaporation'],
        ['constant/reactiveProperties', 'evCoef2', str(evCoef2), 'evaporation'],
        ['constant/reactiveProperties', 'R0', str(R0), 'fermentation'],
        ['constant/reactiveProperties', 'Tm', str(Tm), 'fermentation'],
        ['constant/reactiveProperties', 'deltaT', str(deltaT), 'fermentation'],
        ['constant/reactiveProperties', 'nCoef', str(n), 'evaporation'],
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
        # ['constant/mechanicalProperties', 'nu', str(nu), 'bread'],
        ['constant/mechanicalProperties', 'E', str(E), 'bread'],
        ['constant/mechanicalProperties', 'mu0Raw', str(mu0Raw), 'bread'],
        ['constant/mechanicalProperties', 'kappa0Raw', str(kappa0Raw), 'bread'],
        ['constant/mechanicalProperties', 'muV1Raw', str(muV1Raw), 'bread'],
        ['constant/mechanicalProperties', 'muV2Raw', str(muV2Raw), 'bread'],
        ['constant/mechanicalProperties', 'kappaVRaw', str(kappaV1Raw), 'bread'],
        ['constant/mechanicalProperties', 'mu0Baked', str(mu0Baked), 'bread'],
        ['constant/mechanicalProperties', 'kappa0Baked', str(kappa0Baked), 'bread'],
        ['constant/mechanicalProperties', 'tau1', str(tau1), 'bread'],
        ['constant/mechanicalProperties', 'tau2', str(tau2), 'bread'],
        ['constant/mechanicalProperties', 'tGelat', str(tGelat), 'bread'],
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

if withKynuti:
    baseCase.setParameters(
        [
            ['system/controlDict', 'endTime', str(timeKynuti), ''],
            ['system/controlDict', 'deltaT', '%.5g'%timeStepKynuti, ''],
            ['system/decomposeParDict', 'numberOfSubdomains', str(nCores), ''],
        ]
    )
    baseCase.runCommands(
        [
            'decomposePar > log.decomposePar',
            'foamJob -parallel -screen %s > log.Kynuti' %(dynSolver),
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
        if withKynuti:
            baseCase.setParameters(
                [
                    ['system/controlDict', 'endTime', str(timeKynuti + plusTime1), ''],
                    ['system/controlDict', 'deltaT', '%.5g'%timeStep, ''],
                ]
            )
        else:
            baseCase.runCommands(
                [
                    'decomposePar > log.decomposePar',
                ]
            )
        baseCase.runCommands(
            [
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
                ['system/controlDict', 'endTime', str(timeKynuti + plusTime1 + plusTime2), ''],
                ['constant/transportProperties', 'withDeformation', '0', ''],
                ['system/fvSolution', 'T', str(TRelaxAfter), 'fields'],
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
                'postProcess -func "patchIntegrate(CO2Flux,name=sides)" > log.CO2Flux',
                'rm -rf 0',
                'intMoisture > log.intMoisture',
            ]
        )
    else:
        baseCase.updateTimesParallel()
        baseCase.runCommands(
            [
                'foamJob -parallel -screen postProcess -func "probeZhang" -dict system/probeZhang > log.postProcess',
                'foamJob -parallel -screen postProcess -func "patchIntegrate(CO2Flux,name=sides)" > log.CO2Flux',
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
        
    # -- Load temperature and pressure profiles in probe points
    probesT = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/T'%latestTime, skiprows=3)
    probespG = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/pG'%latestTime, skiprows=3)
    CO2Out = np.loadtxt(baseCase.dir + '/postProcessing/patchIntegrate(CO2Flux,name=sides)/%d/surfaceFieldValue.dat'%latestTime, skiprows=6)

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
    axs[0].plot(probesT[:,0] / 60-timeKynuti/60, probesT[:,1] - 273, 'r', label='center temperature simulation')
    axs[0].plot(probesT[:,0] / 60-timeKynuti/60, probesT[:,2] - 273, 'b', label='center temperature simulation')
    axs[0].set_xlabel("time (min)")
    axs[0].set_ylabel("T (°C)")
    # axs[0].set_xlim(0, 15)
    axs[0].set_title("Temperature evolution in the center and at the surface")
    axs[0].legend()

    # -- Moisture
    axs[1].plot(moistureSim[:,0] / 60 -timeKynuti/60, moistureSim[:,1], 'b', label='simulation')
    axs[1].plot(moistureExp[:,0] , moistureExp[:,1], 'xb', label='experiment')
    axs[1].set_xlabel("time (min)")
    axs[1].set_ylabel("total moisture content (-)")
    # axs[1].set_xlim(0,15)
    axs[1].set_title("Total moisture content in the the bread")
    axs[1].legend()

    axs[2].plot(probespG[:,0] / 60-timeKynuti/60, probespG[:,1], 'r', label='center pressure simulation')
    axs[2].plot(probespG[:,0] / 60-timeKynuti/60, probespG[:,2], 'b', label='surface pressure simulation')
    axs[2].set_xlabel("time (min)")
    axs[2].set_ylabel("p (Pa)")
    # axs[2].set_xlim(0, 15)
    axs[2].set_title("Pressure evolution in the center of the loaf")
    axs[2].legend()

    # -- Displacement
    axs[3].plot(probesT[1:,0] / 60-timeKynuti/60, D[:, 1, 1], 'b', label='simulation DY')
    axs[3].plot(probesT[1:,0] / 60-timeKynuti/60, D[:, 2, 0], 'r', label='simulation DX')
    axs[3].plot(DExp[:,0] / 60, DExp[:,1], 'xr', label='experimental DX')
    axs[3].plot(DExp[:,0] / 60, DExp[:,2], 'xb', label='experimental DY')
    axs[3].set_xlabel("time (min)")
    axs[3].set_ylabel("displacement in X and Y directions")
    # axs[3].set_xlim(0,15)
    axs[3].set_title("Displecement of the bread in vertical (X) and horizontal (Y) directions")
    axs[3].legend()

    # axs[4].plot(CO2Out[:,0] / 60, CO2Out[:,1], 'b', label='CO2 flux simulation')
    # fig.tight_layout()

    # print("CO2Out: ", np.sum(CO2Out[:,1]) * writeInt)

    plt.savefig(baseCase.dir + 'postProcessingPlot.png')
    np.savetxt(os.path.join(baseCase.dir, 'temperatureProbe.dat'), probesT, header='time\tcenter\tsurface\ttop', comments='')
    np.savetxt(os.path.join(baseCase.dir, 'moisture.dat'), moistureSim, header='time\tmoisture', comments='')
                                        
