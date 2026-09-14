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
from compExpSimSingleGraph import saveFigPostProcess

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

nonDeform = True
nonDeform = False

kynuti = False
kynuti = True

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.0e-2 # -- aproximate computational cell size
mSStep = 0.1e-2 # -- aproximate computational cell size

# for expNum in range(len(experiments)): 
for expNum in range(1): 

    '''Internal transport parameters'''
    DFree = 2.6e-5    # -- free volumetric difusivity of the water vapors in CO2 at 300 K
    # Dl = 6e-10  # -- liquid water difusivity in the dough
    # Dl = 1.6e-9  # -- liquid water difusivity in the dough
    # Dl = 5e-11  # -- liquid water difusivity in the dough
    # Dl = 5e-11  # -- liquid water difusivity in the dough
    Dl = 1e-10  # -- liquid water difusivity in the dough
    tortOpen = 2.4   # -- tortuosity
    tortClosed = 70   # -- tortuosity (not used)

    # -- heat conductivity of the dough material with porosity 0, i.e. the 
    # -- absolute term in equation (5) in 
    # -- https://doi.org/10.1016/j.fbp.2008.04.002
    lambdaS = 0.42  # -- heat conductivity of the solid phase (works with addiditional)

    # -- closed-cell bread intristic permeability
    # perm = 1e-15  # -- bread permeability 
    perm = 0.5e-15  # -- bread permeability 
    perm = 0.9e-15  # -- bread permeability 
    # perm = 1e-14  # -- bread permeability 

    # -- heat capacities for the individual phases
    CpS = 1130   # -- solid phase
    CpG = 853  # -- CO2
    CpVapor = 1878 # -- water vapors
    CpL = 4200  # -- liquid phase

    # -- mass density for the individual phases
    rhoS = 701
    # rhoS = 764
    # rhoS = 500
    rhoS = 940
    # rhoS = 970
    # rhoS = 1025

    alphaD0 = 0.84
    if not kynuti:
        alphaD0 = 0.43

    '''Evaporation and CO2 generation parameters'''
    # -- evaporation / condensation coeficient in Hertz-Knudsen equation
    kMPCOpen = 0.01
    kMPCClosed = 0.01

    # -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S) (legacy -- not used)
    evCoef1 = -0.0071
    evCoef2 = 4.5
    n = 0.38

    # -- pre-exponential factor and Tm in CO2 generation kinetics in equation (32) in https://doi.org/10.1002/aic.10518  in (kg/m3/s)
    R0 = 1.8e-3  * 1
    R0 = 2.2e-3  * 1
    # R0 = 1.4e-3  * 1
    # R0 = 0.9e-3
    # R0 = 3e-3
    Tm = 313
    deltaT = 14

    '''Mechanical properties'''
    withDeformation = 1 # -- turn on (1) /off (0) deformation
    if nonDeform:
        withDeformation = 0
    nu = 0.14   # -- Poisson ratio
    # nu = 0.49   # -- Poisson ratio
    E = 30000   # -- Youngs modulus
    # mu0Raw = 230   # kappa = 2*mu*nu/(1-2*nu)  
    mu0Raw = 170   # kappa = 2*mu*nu/(1-2*nu)  
    # mu0Raw = 130   # kappa = 2*mu*nu/(1-2*nu)  
    muV1Raw = 7400 
    # muV1Raw = 4.87e3 
    bakedCoeff = 25
    tau1 = 1.6
    # tau1 = 0.2

    # muV1Raw = 7000
    # tau1 = 1

    '''Numerics and time control'''
    if kynuti:
        timeKynuti = 2400
        # timeKynuti = 200
    else:
        timeKynuti = 200
    timeStepKynuti = 20 # -- computational time step
    timeStepSim = 0.5  # -- computational time step
    timeStepSimNonDef = 0.5  # -- computational time step
    plusTime1 = 450 # -- how long to run with deformation
    # plusTime1 = 800 # -- how long to run with deformation
    plusTime2 = 800 # -- how long to run without deformation

    if nonDeform:
        timeKynuti = 300
        plusTime1 = 1500
        plusTime2 = 0

    writeInt = 30   # -- how often to write results
    writeIntKynuti = 200    # -- how often to write results during kynuti
    nIterKynuti = 200  # -- number of iterations in each time step
    nIterSim = 300  # -- number of iterations in each time step
    nIterSimNonDef = 50  # -- number of iterations in each time step
    dynSolver = 'breadBakingFoam'   # -- used solver
    # dynSolver = 'breadBakingFoamScratch'   # -- used solver
    nCores = 8 # -- number of cores to run the simulation

    # -- relaxation factors
    DRelaxKyn = 0.1
    DFinalRelax = 1

    # -- kynuti
    omegaVRelaxKyn = 1
    omegaCRelaxKyn = 1
    pGRelaxKyn = 1

    # -- deformation simulation
    omegaVRelax = 0.2
    omegaCRelax = 0.2
    pGRelax = 0.2
    DRelax = 1

    # -- non-deformation simulation
    pGRelaxNonDef = 1
    omegaVRelaxNonDef = 0.2
    omegaCRelaxNonDef = 0.2
    TNonDef = 0.1

    '''Boundary and initial conditions'''
    TKynuti = 301
    TStart = 297
    TTop = 200
    TBottom = 200


    kMSidesOmega = 0.01 # -- legacy (not used)
    kMBottomOmega = 0.01    # -- legacy (not used)
    kMTop = 3e-3   # -- legacy (not used) external mass transfer coeficient
    alphaG = 18 # -- external heat transfer coeficient 
    alphaGBottom = 18 # -- external heat transfer coeficient 

    alphaG = 10 # -- external heat transfer coeficient 
    alphaGBottom = 10
    
     # -- external heat transfer coeficient 

    '''Post-processing'''
    fig, axs = plt.subplots(1, 1, figsize=(16, 9))  # figure with plots

    outFolder = '../ZZ_cases/2026/V30/exp%d_nonDef_%s/V19_mechFrompaper_mu0_%g_Dl_%g_kOp_%g_kCl_%g_nu_%g_mu0_%g_muV1_%g_mS_%g_lambda_%g_kH_%g_kHB_%g_r0_%g_per_%g/' % (expNum, str(nonDeform), mu0Raw, Dl, kMPCOpen, kMPCClosed, nu, mu0Raw, muV1Raw, mSStep, lambdaS, alphaG, alphaGBottom, R0, perm)
    # baseCaseDir = '../ZZ_cases/2026/V30/exp%d_nonDef_%s/V23_losingMoisture_Dl_%g_kOp_%g_kCl_%g_nu_%g_mu0_%g_muV1_%g_mS_%g_lambda_%g_kH_%g_kHB_%g_r0_%g_per_%g/' % (expNum, str(nonDeform), Dl, kMPCOpen, kMPCClosed, nu, mu0Raw, muV1Raw, mSStep, lambdaS, alphaG, alphaGBottom, R0, perm)


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
    tau0 = 10

    # -- prepare blockMeshDict using luckas python class
    if prepBlockMesh:   
        # prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=experiments[expNum]['up'])
        if kynuti:
            prep3DMeshOurExp(experiments[expNum]['rLoafKynuti'], experiments[expNum]['rLoafKynuti'], experiments[expNum]['hLoafKynuti'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=5e-3)
        else:
            prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=5e-3)

        if nonDeform:
            # prep2DMeshOurExp(experiments[expNum]['rLoaf']*2, experiments[expNum]['hLoaf']*2, x0, y0, z0, dA, dX, dY, dZ, grX, grY, grZ, baseCase, for3D=True)
            prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf']+2.5e-2, dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=0e-2,nonDeform=nonDeform)
        # prep3DMeshOurExp(experiments[expNum]['rLoafKynuti'], experiments[expNum]['rLoafKynuti'], experiments[expNum]['hLoafKynuti'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=False, up=0)
        # prep3DMeshForSnappy(1e-4, rLoaf1, 1e-4, rLoaf2, 1e-4, hLoaf, dA, baseCase, for2DExtrude=False)

    # CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
    # 1) BOUNDARY CONDITIONS
    # -- change in tutorial case
    baseCase.setParameters(
        [
            ['0.org/omegaV', 'kM', str(kMSidesOmega), 'sides'],
            ['0.org/omegaV', 'kM', str(kMBottomOmega), 'bottom'],
            ['0.org/omegaV', 'kM', str(kMBottomOmega), 'bottom2'],
            ['0.org/omegaC', 'kM', str(kMSidesOmega), 'sides'],
            ['0.org/omegaC', 'kM', str(kMBottomOmega), 'bottom'],
            ['0.org/omegaC', 'kM', str(kMBottomOmega), 'bottom2'],
            # ['0.org/pG', 'kM', str(kMSides), 'sides'],
            # ['0.org/pG', 'kM', str(kMBottom), 'bottom'],
            # ['0.org/pG', 'kM', str(kMBottom), 'bottom2'],
            ['0.org/T', 'alpha', str(alphaG), 'sides'],
            ['0.org/T', 'alpha', str(alphaGBottom), 'bottom'],
            ['0.org/T', 'alpha', str(alphaGBottom), 'bottom2'],
            # ['0.org/alphaL', 'internalField', 'uniform %g'% (alphaL), ''],
            # ['0.org/alphaS', 'internalField', 'uniform %g'% (alphaS), ''],
            ['0.org/T', 'internalField', 'uniform %g'% (TStart), ''],
        ]
    )

    # -- get the external temperature
    cleanDataInRange = loadDataFrame(expDir, experiments[expNum])
    bakingCurve = getTemps(cleanDataInRange, 4)

    # -- change external temperature
    with open(os.path.join(baseCase.dir, "constant", "TInfTable"), "w") as fl:
        fl.writelines("(\n")
        # if not nonDeform:
        fl.writelines("\t(0\t%f)\n"%TKynuti)
        fl.writelines("\t(%d\t%f)\n"%(timeKynuti, TKynuti))
        bakingCurve[:, 1] = TTop
        for i in range(bakingCurve.shape[0]):
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeKynuti+0.1, bakingCurve[i,1]))
            fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeKynuti+0.1, bakingCurve[i,1] + 273.15))
        # else:
        #     for i in range(bakingCurve.shape[0]):
        #         fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60, bakingCurve[i,1] + 273.15))
        fl.writelines(")\n")

    with open(os.path.join(baseCase.dir, "constant", "TInfTableBottom"), "w") as fl:
        fl.writelines("(\n")
        # if not nonDeform:
        fl.writelines("\t(0\t%f)\n"%TKynuti)
        fl.writelines("\t(%d\t%f)\n"%(timeKynuti, TKynuti))
        bakingCurve[:, 1] = TBottom
        for i in range(bakingCurve.shape[0]):
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeKynuti+0.1, bakingCurve[i,1]))
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeKynuti+0.1, 10 + bakingCurve[i,1] + 273.15))
            fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeKynuti+0.1, bakingCurve[i,1] + 273.15))
        # else:
        #     for i in range(bakingCurve.shape[0]):
        #         fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60, bakingCurve[i,1] + 273.15))
        fl.writelines(")\n")

    # 2) constant/transportProperties
    baseCase.setParameters(
        [
            ['constant/transportProperties', 'withDeformation', str(withDeformation), ''],
            ['constant/transportProperties', 'permGLViscG', str(perm), ''],
            ['constant/transportProperties', 'tortOpen', str(tortOpen), ''],
            ['constant/transportProperties', 'tortClosed', str(tortClosed), ''],
            ['constant/transportProperties', 'alphaD0', str(alphaD0), ''],
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
            # ['constant/reactiveProperties', 'TGelEv', str(tGelatEv), 'gelatization'],
            # ['constant/reactiveProperties', 'alphaKept', str(alphaKept), 'gelatization']
        ]
    )
            
    # 5 system/controlDict
    baseCase.setParameters(
        [
            # ['system/controlDict', 'endTime', str(plusTime1), ''],
            ['system/controlDict', 'endTime', str(timeKynuti), ''],
            ['system/controlDict', 'deltaT', '%.5g'%timeStepKynuti, ''],
            ['system/controlDict', 'writeInterval', '%.5g'%writeIntKynuti, ''],
        ]
    )

    # 6) fvSolutions
    baseCase.setParameters(
        [
            ['system/fvSolution', 'nOuterCorrectors', str(nIterKynuti), 'PIMPLE'],
            ['system/fvSolution', 'D', str(DRelaxKyn), 'fields'],
            ['system/fvSolution', 'DFinal', str(DFinalRelax), 'fields'],
            ['system/fvSolution', 'omegaV', str(omegaVRelaxKyn), 'fields'],
            ['system/fvSolution', 'omegaC', str(omegaCRelaxKyn), 'fields'],
            ['system/fvSolution', 'pG', str(pGRelaxKyn), 'fields'],
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
            ]
        )
        # if nonDeform:
        #     baseCase.runCommands(
        #         [
        #             'snappyHexMesh -overwrite > log.snappyHexMesh',
        #         ]
        #     )
        baseCase.runCommands(
            [   
                # 'snappyHexMesh -overwrite > log.snappyHexMesh',
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

        if plusTime1 > 0:
            baseCase.setParameters(
                [
                    ['system/controlDict', 'endTime', str(timeKynuti + plusTime1), ''],
                    ['system/controlDict', 'deltaT', '%.5g'%timeStepSim, ''],
                    # ['constant/transportProperties', 'withDeformation', '0', ''],
                    ['system/fvSolution', 'nOuterCorrectors', str(nIterSim), 'PIMPLE'],
                    ['system/fvSolution', 'omegaV', str(omegaVRelax), 'fields'],
                    ['system/fvSolution', 'omegaC', str(omegaCRelax), 'fields'],
                    ['system/fvSolution', 'D', str(DRelax), 'fields'],
                    ['system/fvSolution', 'pG', str(pGRelax), 'fields'],
                    ['system/controlDict', 'writeInterval', '%.5g'%writeInt, ''],


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

        # -- run the rest of the simualation without further deformation
        if plusTime2 > 0:
            baseCase.setParameters(
                [
                    ['system/controlDict', 'endTime', str(timeKynuti + plusTime1 + plusTime2), ''],
                    ['system/controlDict', 'deltaT', '%.5g'%timeStepSimNonDef, ''],
                    ['constant/transportProperties', 'withDeformation', '0', ''],
                    ['system/fvSolution', 'omegaV', str(omegaVRelaxNonDef), 'fields'],
                    ['system/fvSolution', 'omegaC', str(omegaCRelaxNonDef), 'fields'],
                    ['system/fvSolution', 'T', str(TNonDef), 'fields'],
                    ['system/fvSolution', 'pG', str(pGRelaxNonDef), 'fields'],
                    ['system/fvSolution', 'nOuterCorrectors', str(nIterSimNonDef), 'PIMPLE'],

                ]
            )
            if nCores > 1:
                baseCase.runCommands(
                    [
                        'foamJob -parallel -screen %s > log.%s_3' %(dynSolver,dynSolver),
                    ]
                )
            else:
                baseCase.runCommands(
                    [
                        '%s > log.%s_3' %(dynSolver,dynSolver),
                    ]
                )


    # # POST-PROCESSING=======================================================
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
                    'getBoundPoints > log.getBoundPoints'
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
                    'foamJob -parallel -screen intWeigth > log.intWeigth',
                    'foamJob -parallel -screen getBoundPoints > log.getBoundPoints'
                ]
            )
            for i in range(len(experiments[expNum]['probes'])):
                # thermoOffsetCorr = np.array([experiments[expNum]['thermoOffset'][0], experiments[expNum]['thermoOffset'][2], -experiments[expNum]['thermoOffset'][1]])
                thermoOffsetCorr = experiments[expNum]['thermoOffset']
                probesCorr = experiments[expNum]['probes'][i] - thermoOffsetCorr
                probesCorr[1] = np.sqrt(probesCorr[1]**2 + probesCorr[2]**2)
                probesCorr[2] = 0
                baseCase.runCommands(
                    [
                        'foamJob -parallel -screen TLFProbe -point "(%.5g %.5g %.5g)" > log.TPoint%d' %(probesCorr[0], probesCorr[1], probesCorr[2], i+1),
                    ]
                )

        if nCores > 1:
            latestTime = baseCase.latestParTime
        else:
            latestTime  = baseCase.latestTime

        shapeBef = np.loadtxt(baseCase.dir + 'Shape/D_values_%d.dat' % timeKynuti, skiprows=1)   
        shapeAft = np.loadtxt(baseCase.dir + 'Shape/D_values_%d.dat' % latestTime, skiprows=1)   
        shapeBef = shapeBef[np.argsort(shapeBef[:, 0])]
        shapeAft = shapeAft[np.argsort(shapeAft[:, 0])]
        np.savetxt(os.path.join(baseCase.dir, 'Shape/D_values_%d_sorted.dat' % timeKynuti), shapeBef, header='x\ty\tz', comments='')
        np.savetxt(os.path.join(baseCase.dir, 'Shape/D_values_%d_sorted.dat' % latestTime), shapeAft, header='x\ty\tz', comments='')

        # -- gather the displacement data from probe points
        # rows = []
        # lines = []
        # D = []
        # nProbes = 2
        # if nCores > 1:
        #     latestTime = baseCase.latestParTime
        # else:
        #     latestTime  = baseCase.latestTime
        # with open(baseCase.dir + '/postProcessing/probeOur/%d/D'%latestTime, 'r') as fl:
        #     lines = fl.readlines()
        #     lines = lines[nProbes+1:]
        #     # print(lines)
        #     for line in lines:
        #         parts = line.split(") (")
        #         first_entry = parts[0].split(maxsplit=1)
        #         vectors = [first_entry[1]] if len(first_entry) > 1 else []
        #         vectors.extend(parts[1:])

        #         vectors = [
        #             tuple(map(float, vec.replace("(", "").replace(")", "").split()))
        #             for vec in vectors
        #         ]
        #         rows.append(vectors)

        # # -- Convert displacements to numpy array
        # D = np.array(rows)
        # np.savetxt(os.path.join(baseCase.dir, 'sim_DX_DY.dat'), np.column_stack([D[:,0, 0], D[:,1,1]]), header='DX DY', comments='')
            
        # # -- Load temperature profiles in probe points
        # # probesT = np.loadtxt(baseCase.dir + '/postProcessing/probeZhang/%d/T'%latestTime, skiprows=3)

        # # -- Load total moisture evolution 
        # TPoint6 = readDataFromLogFile("%s/log.TPoint6" %baseCase.dir)
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
        # axs.plot(TPoint6[:,0] / 60 , D[:, 0, 0], 'b', label='simulation DX')
        # axs.plot(TPoint6[:,0] / 60 , D[:, 1, 1], 'r', label='simulation DY')
        # axs.plot(experiments[expNum]['expDispl'][:,0], experiments[expNum]['expDispl'][:,1], 'xb', label='experimental DX')
        # axs.plot(experiments[expNum]['expDispl'][:,0], experiments[expNum]['expDispl'][:,2], 'xr', label='experimental DY')
        # # axs[2].plot(DExp[:,0] / 60, DExp[:,2], 'xb', label='experimental DY')
        # axs.set_xlabel("time (min)")
        # axs.set_ylabel("displacement in X and Y directions")
        # # axs[2,0].set_xlim(0,35)
        # axs.set_title("Displacement of the bread in vertical (X) and horizontal (Y) directions")
        # axs.legend()
        # fig.tight_layout()

        # plt.savefig(baseCase.dir + 'postProcessingPlot.png')
                                            
        saveFigPostProcess(timeKynuti, outFolder)
