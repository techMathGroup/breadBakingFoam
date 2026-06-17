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
from expDict import *
import os

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/bread3DOurExp/' # -- base case for simulation
outFolder = '../ZZ_cases/2026/test/'
expDir = os.path.join('..', 'Experiments2026')

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
# prepBlockMesh = False    # -- preparation of the blockMeshDict script
# makeGeom = False # -- creation of the geometry for computation
# runDynSim = False    # -- run simulation
runPostProcess = True   # -- run post-processing

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.08e-2 # -- aproximate computational cell size
# rLoaf1 = 8.5e-2  # -- loaf radius                
# rLoaf2 = 8.5e-2  # -- loaf radius                
# hLoaf = 7e-2  # -- loaf height 
# up = 1e-2

expNum = 0

'''Internal transport parameters'''
# -- free volumetric difusivity of the water vapors in CO2 at 300 K
DFree = 5e-6 

# -- heat conductivity of the dough material with porosity 0, i.e. the 
# -- absolute term in equation (5) in 
# -- https://doi.org/10.1016/j.fbp.2008.04.002
lambdaS = 0.55

perm = 2.6e-12  # -- bread permeability 

# -- heat capacities for the individual phases
CpS = 1450   # -- solid phase
CpG = 853  # -- CO2
CpVapor = 1878 # -- water vapors
CpL = 4200  # -- liquid phase

# -- mass densities for the individual phases
rhoS = 1200  # -- solid density    
rhoL = 1000  # -- liquid density   
alphaL = 0.24
alphaS = 0.28

'''Evaporation and CO2 generation parameters'''
# -- evaporation / condensation coeficient in Hertz-Knudsen equation
kMPC = 0.04

# -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S)
evCoef1 = -0.0071
evCoef2 = 4.5
n = 0.38
alphaKept = 0.25

# -- pre-exponential factor and Tm in CO2 generation kinetics 
# -- in equation (32) in https://doi.org/10.1002/aic.10518
R0 = 20e-4 
# R0 = 0 
Tm = 310
tGelatEv = 57
Tm = 310

'''Mechanical properties'''
withDeformation = 1 # -- turn on (1) /off (0) deformation
nu = 0.15   # -- Poisson ratio
E = 12000   # -- Youngs modulus
tau0 = 2

'''Numerics'''
timeStep = 1  # -- computational time step
plusTime1 = 1800 # -- how long to run with deformation
plusTime2 = 0 # -- how long to run without deformation
writeInt = 30   # -- how often to write results
nIter = 30  # -- number of iterations in each time step
dynSolver = 'breadBakingFoam'   # -- used solver
nCores = 8 # -- number of cores to run the simulation

# -- relaxation factors
DRelax = 0.1
DFinalRelax = 1

'''Boundary conditions'''
kMSides = 1e-3   # -- external mass transfer coeficient
kMBottom = 0.01   # -- external mass transfer coeficient
kMBottom2 = 1e-3
kMTop = 3e-3   # -- external mass transfer coeficient
alphaG = 20 # -- external heat transfer coeficient 

'''Post-processing'''
fig, axs = plt.subplots(1, 1, figsize=(16, 9))  # figure with plots

# baseCaseDir = '../ZZ_cases/02_ourBig/bread3DOurExp_newUpdates_Cps1450_alphaG_%g_kmsides%g_kmbottom%g_kmtop%g_r0_%g_perm_%g/' % (alphaG, kMSides, kMBottom, kMTop, R0, perm)
# baseCaseDir = '../ZZ_cases/2026/exp%d/alphaG_%g_kmsides%g_kmbottom%g_kmtop%g_r0_%g_perm_%g/' % (expNum, alphaG, kMSides, kMBottom, kMTop, R0, perm)
outFolder = '../ZZ_cases/2026/exp%d_newTest_newBound/viceVersa_tauLow_pG_alphaG_%g_kmsides%g_kmbottom%g_kmtop%g_r0_%g_perm_%g/' % (expNum, alphaG, kMSides, kMBottom, kMTop, R0, perm)


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
    prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=experiments[expNum]['up'])
    # prep3DMeshForSnappy(1e-4, rLoaf1, 1e-4, rLoaf2, 1e-4, hLoaf, dA, baseCase, for2DExtrude=False)

# CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
# 1) BOUNDARY CONDITIONS
# -- change in tutorial case
baseCase.setParameters(
    [
        ['0.org/omegaV', 'kM', str(kMSides), 'sides'],
        ['0.org/omegaV', 'kM', str(kMBottom), 'bottom'],
        ['0.org/omegaV', 'kM', str(kMBottom2), 'bottom2'],
        ['0.org/pG', 'kM', str(kMSides), 'sides'],
        ['0.org/pG', 'kM', str(kMBottom), 'bottom'],
        ['0.org/pG', 'kM', str(kMBottom2), 'bottom2'],
        ['0.org/T', 'alpha', str(alphaG), 'sides'],
        ['0.org/T', 'alpha', str(alphaG), 'bottom'],
        ['0.org/T', 'alpha', str(alphaG), 'bottom2'],
        ['0.org/alphaL', 'internalField', 'uniform %g'% (alphaL), ''],
        ['0.org/alphaS', 'internalField', 'uniform %g'% (alphaS), ''],
    ]
)

# -- get the external temperature
cleanDataInRange = loadDataFrame(expDir, experiments[expNum])
bakingCurve = getTemps(cleanDataInRange, 4)


# -- change external temperature
with open(os.path.join(baseCase.dir, "constant", "TInfTable"), "w") as fl:
    fl.writelines("(\n")
    for i in range(bakingCurve.shape[0]):
        fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60, bakingCurve[i,1] + 273.15))
    fl.writelines(")\n")

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
        ['constant/reactiveProperties', 'nCoef', str(n), 'evaporation'],
        ['constant/reactiveProperties', 'TGelEv', str(tGelatEv), 'gelatization'],
        ['constant/reactiveProperties', 'alphaKept', str(alphaKept), 'gelatization']
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
baseCase.addToDictionary(
    [
        ['constant/mechanicalProperties','tau0 %g;\n'%tau0, ''],
    ]
)

# -- prepare geom
if makeGeom:
    baseCase.runCommands(
        [
            'chmod 755 ./* -R',
            'blockMesh > log.blockMesh',
            'extrudeMesh > log.extrudeMesh',
            # 'topoSet > log.topoSet',
            # 'createPatch -overwrite > log.createPatch',
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
    expDataDispl = np.loadtxt(baseCaseDir + 'ZZ_dataForPostProcessing/exp_DX_DY.dat', skiprows=1)
    
    # -- run post-processing tasks
    if nCores == 1:
        baseCase.updateTimes()
        baseCase.runCommands(
            [
                'postProcess -func "probeOur" -dict system/probeOur > log.postProcess',
                'TLFProbe -point "(0.012 1e-4 1e-4)" > log.TPoint6',
                'TLFProbe -point "(0.061 1e-3 1e-3)" > log.TPoint7',
                'TLFProbe -point "(0.027 0.047 1e-4)" > log.TPoint5',
                'TLFProbe -point "(0.032 0.041 1e-4)" > log.TPoint8',
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
                'foamJob -parallel -screen TLFProbe -point "(1e-3 1e-3 0)" > log.TPoint6',
                # 'foamJob -parallel -screen TLFProbe -point "(0.061 1e-3 0)" > log.TPoint7',
                # 'foamJob -parallel -screen TLFProbe -point "(0.027 0.047 0)" > log.TPoint5',
                # 'foamJob -parallel -screen TLFProbe -point "(0.032 0.041 0)" > log.TPoint8',
                # 'foamJob -parallel -screen TLFProbe -point "(0.022 1e-4 0)" > log.TPoint66',
                # 'foamJob -parallel -screen TLFProbe -point "(0.071 1e-3 0)" > log.TPoint77',
                # 'foamJob -parallel -screen TLFProbe -point "(0.037 0.047 0)" > log.TPoint55',
                # 'foamJob -parallel -screen TLFProbe -point "(0.042 0.041 0)" > log.TPoint88',
                'foamJob -parallel -screen intMoisture > log.intMoisture',
            ]
        )
        for i in range(len(experiments[expNum]['probes'])):
            baseCase.runCommands(
                [
                    'foamJob -parallel -screen TLFProbe -point "(%.5g %.5g %.5g)" > log.TPoint%d' %(experiments[expNum]['probes'][i][0], experiments[expNum]['probes'][i][1], experiments[expNum]['probes'][i][2], i+1),
                ]
            )

    # -- gather the displacement data from probe points
    rows = []
    lines = []
    D = []
    nProbes = 2
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
    np.savetxt(os.path.join(baseCase.dir, 'sim_DX_DY.dat'), np.column_stack([D[:,0, 0], D[:,1,1]]), header='DX DY', comments='')
        
    # -- Load temperature profiles in probe points
    # probesT = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/T'%latestTime, skiprows=3)

    # -- Load total moisture evolution 
    TPoint6 = readDataFromLogFile("%s/log.TPoint6" %baseCase.dir)
    # TPoint7 = readDataFromLogFile("%s/log.TPoint7" %baseCase.dir)
    # TPoint5 = readDataFromLogFile("%s/log.TPoint5" %baseCase.dir)
    # TPoint8 = readDataFromLogFile("%s/log.TPoint8" %baseCase.dir)

    # TPoint66 = readDataFromLogFile("%s/log.TPoint66" %baseCase.dir)
    # TPoint77 = readDataFromLogFile("%s/log.TPoint77" %baseCase.dir)
    # TPoint55 = readDataFromLogFile("%s/log.TPoint55" %baseCase.dir)
    # TPoint88 = readDataFromLogFile("%s/log.TPoint88" %baseCase.dir)
    # moistureSim = readDataFromLogFile("%s/log.intMoisture" %baseCase.dir)


    # -- Temperatures
    # axs[0,0].plot(expData[:,-2],expData[:,0], '--r',  label='exp. point 5')
    # axs[0,0].plot(expData[:,-2],expData[:,1], '--g',  label='exp. point 6')
    # axs[0,0].plot(expData[:,-2],expData[:,2], '--b',  label='exp. point 7')
    # axs[0,0].plot(expData[:,-2],expData[:,3], '--m',  label='exp. point 8')
    # axs[0,1].plot(expData2[:,-2],expData2[:,0], '--r',  label='exp. point 5')
    # axs[0,1].plot(expData2[:,-2],expData2[:,1], '--g',  label='exp. point 6')
    # axs[0,1].plot(expData2[:,-2],expData2[:,2], '--b',  label='exp. point 7')
    # axs[0,1].plot(expData2[:,-2],expData2[:,3], '--m',  label='exp. point 8')
    # axs[0,0].plot(expData[:,-2],expData[:,0], '--r',  label='exp. point 5')
    # axs[0,0].plot(expData[:,-2],expData[:,1], '--g',  label='exp. point 6')
    # axs[0,0].plot(expData[:,-2],expData[:,2], '--b',  label='exp. point 7')
    # axs[0,0].plot(expData[:,-2],expData[:,3], '--m',  label='exp. point 8')
    # axs[0,1].plot(expData2[:,-2],expData2[:,0], '--r',  label='exp. point 5')
    # axs[0,1].plot(expData2[:,-2],expData2[:,1], '--g',  label='exp. point 6')
    # # axs[0,1].plot(expData2[:,-2],expData2[:,2], '--b',  label='exp. point 7')
    # axs[0,1].plot(expData2[:,-2],expData2[:,3], '--m',  label='exp. point 8')

    # axs[0].plot(TExpPoint2[:,-1],TExpPoint2[:,0], '--b',  label='exp. point 2')
    # axs[0].plot(TExpPoint3[:,-1],TExpPoint3[:,0], '--g',  label='exp. point 4')
    # axs[0].plot(TExpSurface[:,0],TExpSurface[:,1], 'xb', label='surface temperature experiment')
    # axs[0,0].plot(TPoint5[:,0] / 60, TPoint5[:,1] - 273, 'r', label='sim. point 5')
    # axs[0,0].plot(TPoint6[:,0] / 60, TPoint6[:,1] - 273, 'g', label='sim. point 6')
    # axs[0,0].plot(TPoint7[:,0] / 60, TPoint7[:,1] - 273, 'b', label='sim. point 7')
    # axs[0,0].plot(TPoint8[:,0] / 60, TPoint8[:,1] - 273, 'm', label='sim. point 8')

    # axs[0,1].plot(TPoint55[:,0] / 60 , TPoint55[:,1] - 273, 'r', label='sim. point 5')
    # axs[0,1].plot(TPoint66[:,0] / 60 , TPoint66[:,1] - 273, 'g', label='sim. point 6')
    # axs[0,1].plot(TPoint77[:,0] / 60 , TPoint77[:,1] - 273, 'b', label='sim. point 7')
    # axs[0,1].plot(TPoint88[:,0] / 60 , TPoint88[:,1] - 273, 'm', label='sim. point 8')
    # axs[0,1].plot(TPoint5[:,0] / 60, TPoint5[:,1] - 273, 'r', label='sim. point 5')
    # axs[0,1].plot(TPoint6[:,0] / 60, TPoint6[:,1] - 273, 'g', label='sim. point 6')
    # axs[0,1].plot(TPoint7[:,0] / 60, TPoint7[:,1] - 273, 'b', label='sim. point 7')
    # axs[0,1].plot(TPoint8[:,0] / 60, TPoint8[:,1] - 273, 'm', label='sim. point 8')
    # axs[0].plot(probesT[:,0] / 60, probesT[:,2] - 273, 'b', label='center temperature simulation')
    # axs[0,0].set_xlabel("time (min)")
    # axs[0,0].set_ylabel("T (°C)")
    # axs[0,0].set_ylim(20,120)
    # axs[0,0].set_xlim(0, 35)
    # axs[0,0].set_title("Temperature evolution in the center and at the surface")
    # axs[0,0].legend()

    # axs[0,1].set_xlabel("time (min)")
    # axs[0,1].set_ylabel("T (°C)")
    # axs[0,1].set_ylim(20,120)
    # axs[0,1].set_xlim(0, 35)
    # axs[0,1].set_title("Temperature evolution in the center and at the surface")
    # axs[0,1].legend()

    # # -- Moisture
    # axs[1,0].plot(moistureSim[:,0] / 60 , moistureSim[:,1], 'b', label='simulation')
    # axs[1,0].plot(expData[:,-2], expData[:,-1], '--b', label='experiment')
    # axs[1,0].set_xlabel("time (min)")
    # axs[1,0].set_ylabel("total moisture content (-)")
    # axs[1,0].set_ylim(0.35,0.7)
    # # axs[1].set_xlim(0,28)
    # axs[1,0].set_title("Total moisture content in the the bread")
    # axs[1,0].legend()

    # axs[1,1].plot(moistureSim[:,0] / 60 , moistureSim[:,1], 'b', label='simulation')
    # axs[1,1].plot(expData2[:,-2], expData2[:,-1], '--b', label='experiment')
    # axs[1,1].set_xlabel("time (min)")
    # axs[1,1].set_ylabel("total moisture content (-)")
    # axs[1,1].set_ylim(0.35,0.7)
    # # axs[1].set_xlim(0,28)
    # axs[1,1].set_title("Total moisture content in the the bread")
    # axs[1,1].legend()

    # # -- Displacement
    axs.plot(TPoint6[:,0] / 60 , D[:, 0, 0], 'b', label='simulation DX')
    axs.plot(TPoint6[:,0] / 60 , D[:, 1, 1], 'r', label='simulation DY')
    axs.plot(experiments[expNum]['expDispl'][:,0], experiments[expNum]['expDispl'][:,1], 'xb', label='experimental DX')
    axs.plot(experiments[expNum]['expDispl'][:,0], experiments[expNum]['expDispl'][:,2], 'xr', label='experimental DY')
    # axs[2].plot(DExp[:,0] / 60, DExp[:,2], 'xb', label='experimental DY')
    axs.set_xlabel("time (min)")
    axs.set_ylabel("displacement in X and Y directions")
    # axs[2,0].set_xlim(0,35)
    axs.set_title("Displacement of the bread in vertical (X) and horizontal (Y) directions")
    axs.legend()
    fig.tight_layout()

    plt.savefig(baseCase.dir + 'postProcessingPlot.png')
                                        
