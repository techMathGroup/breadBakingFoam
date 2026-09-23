#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python script to set up and run bread baking simulations according to 
# custom experiments

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
from compExpSimSingleGraphFromDat import saveFigPostProcess

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/breadAx2DOurExp/' # -- base case for simulation
outFolder = '../ZZ_cases/01_breadAx2DOurExp3/'
# expDir = os.path.join('..', 'Experiments2026') # -- when comparing experiments

# WHAT SHOULD RUN=======================================================
prepBlockMesh = True    # -- preparation of the blockMeshDict script
makeGeom = True # -- creation of the geometry for computation
runDynSim = True    # -- run simulation
# prepBlockMesh = False    # -- preparation of the blockMeshDict script
# makeGeom = False # -- creation of the geometry for computation
# runDynSim = False    # -- run simulation
runPostProcess = True   # -- run post-processing
runWithSlurm = True

proofing = False  # -- proofing included
proofing = True  # -- proofing included

# DEFINE PARAMETERS=====================================================
'''Geometry parameters'''
mSStep = 0.1e-2 # -- aproximate computational cell size

# for expNum in range(len(experiments)): 
for expNum in range(1): 

    '''Internal transport parameters'''
    DFree = 2.6e-5    # -- free volumetric difusivity of the water vapors in CO2 at 300 K
    Dl = 6e-11  # -- liquid water difusivity in the dough
    tortOpen = 2.4   # -- tortuosity
    tortClosed = 70   # -- tortuosity (not used)

    # -- heat conductivity of the dough material with porosity 0, i.e. the 
    # -- absolute term in equation (5) in 
    # -- https://doi.org/10.1016/j.fbp.2008.04.002
    lambdaS = 0.42  # -- heat conductivity of the solid phase (works with addiditional)

    # -- intrinsic gas permeabilities of raw and baked bread
    gasPermeabilityRaw = 0.9e-15
    gasPermeabilityBaked = 1e-11

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

    '''Evaporation and CO2 generation parameters'''
    # -- evaporation / condensation coeficient in Hertz-Knudsen equation
    kMPCOpen = 0.01
    kMPCClosed = 0.01

    # -- parameters for Oswin model (https://doi.org/10.1016/0260-8774(91)90020-S) (legacy -- not used)
    evCoef1 = -0.0071
    evCoef2 = 4.5
    n = 0.38

    # -- pre-exponential factor and Tm in CO2 generation kinetics in equation (32) in https://doi.org/10.1002/aic.10518  in (kg/m3/s)
    R0 = 2.3e-3  
    # R0 = 1.8e-3  
    Tm = 313
    deltaT = 14

    '''Mechanical properties'''
    withDeformation = 1 # -- turn on (1) /off (0) deformation

    nu = 0.14   # -- Poisson ratio
    E = 30000   # -- Youngs modulus (legacy -- not used) 
    mu0Raw = 147   # kappa = 2*mu*nu/(1-2*nu)   
    muV1Raw = 7400 
    bakedCoeff = 25
    tau1 = 1.6

    '''Numerics and time control'''
    if proofing:
        timeProofing = 2400
    else:
        timeProofing = 200
    timeStepProofing = 20 # -- computational time step for proofing
    timeStepSim = 0.5  # -- computational time step for deformable simulation
    timeStepSimNonDef = 0.5  # -- computational time step non-deformable simulation
    plusTime1 = 450 # -- how long to run with deformation 
    plusTime2 = 750 # -- how long to run without deformation

    writeInt = 30   # -- how often to write results
    writeIntProofing = 200    # -- how often to write results during proofing
    nIterProofing = 300  # -- number of iterations in each time step
    nIterSim = 250  # -- number of iterations in each time step
    nIterSimNonDef = 50  # -- number of iterations in each time step
    dynSolver = 'breadBakingFoam'   # -- used solver
    # dynSolver = 'breadBakingFoamScratch'   # -- used solver
    nCores = 8 # -- number of cores to run the simulation

    # -- relaxation factors
    DRelaxKyn = 0.1
    DFinalRelax = 1

    # -- proofing
    omegaVRelaxKyn = 1
    omegaCRelaxKyn = 1
    pGRelaxKyn = 1

    # -- deformation simulation
    omegaVRelax = 0.2
    omegaCRelax = 0.2
    pGRelax = 0.2
    DRelax = 0.8
    TRelaxKyn = 0.2

    # -- non-deformation simulation
    pGRelaxNonDef = 1
    omegaVRelaxNonDef = 0.2
    omegaCRelaxNonDef = 0.2
    TNonDef = 0.1

    '''Boundary and initial conditions'''
    TProofing = 301
    TStart = 298
    TTop = 200
    TBottom = 200


    kMSidesOmega = 0.01 # -- external mass transfer coeficient 
    kMBottomOmega = 0.01

    alphaG = 10 # -- external heat transfer coeficient 
    alphaGBottom = 14
    
     # -- external heat transfer coeficient 

    '''Post-processing'''
    fig, axs = plt.subplots(1, 1, figsize=(16, 9))  # figure with plots

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
    node = "kraken-x3"

    # -- prepare blockMeshDict using luckas python class
    if prepBlockMesh:   
        # prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=experiments[expNum]['up'])
        if proofing:
            prep3DMeshOurExp(experiments[expNum]['rLoafKynuti'], experiments[expNum]['rLoafKynuti'], experiments[expNum]['hLoafKynuti'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=5e-3)
        else:
            prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=5e-3)

        # if nonDeform:
        #     # prep2DMeshOurExp(experiments[expNum]['rLoaf']*2, experiments[expNum]['hLoaf']*2, x0, y0, z0, dA, dX, dY, dZ, grX, grY, grZ, baseCase, for3D=True)
        #     prep3DMeshOurExp(experiments[expNum]['rLoaf'], experiments[expNum]['rLoaf'], experiments[expNum]['hLoaf']+2.5e-2, dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=True, up=0e-2,nonDeform=nonDeform)
        # prep3DMeshOurExp(experiments[expNum]['rLoafKynuti'], experiments[expNum]['rLoafKynuti'], experiments[expNum]['hLoafKynuti'], dX, dY, dZ, grX, grY, grZ, baseCase, for2DExtrude=False, up=0)
        # prep3DMeshForSnappy(1e-4, rLoaf1, 1e-4, rLoaf2, 1e-4, hLoaf, dA, baseCase, for2DExtrude=False)

    # CHANGE THE PARAMETERS IN OPENFOAM DICTIONARIES========================
    # 1) BOUNDARY CONDITIONS
    # -- change in tutorial case
    baseCase.setParameters(
        [
            ['0.org/omegaV', 'kM', str(kMSidesOmega), 'sides'],
            ['0.org/omegaV', 'kM', str(kMBottomOmega), 'bottom'],
            # ['0.org/omegaV', 'kM', str(kMBottomOmega), 'bottom2'],
            ['0.org/omegaC', 'kM', str(kMSidesOmega), 'sides'],
            ['0.org/omegaC', 'kM', str(kMBottomOmega), 'bottom'],
            # ['0.org/omegaC', 'kM', str(kMBottomOmega), 'bottom2'],
            # ['0.org/pG', 'kM', str(kMSides), 'sides'],
            # ['0.org/pG', 'kM', str(kMBottom), 'bottom'],
            # ['0.org/pG', 'kM', str(kMBottom), 'bottom2'],
            ['0.org/T', 'alpha', str(alphaG), 'sides'],
            ['0.org/T', 'alpha', str(alphaGBottom), 'bottom'],
            # ['0.org/T', 'alpha', str(alphaGBottom), 'bottom2'],
            # ['0.org/alphaL', 'internalField', 'uniform %g'% (alphaL), ''],
            # ['0.org/alphaS', 'internalField', 'uniform %g'% (alphaS), ''],
            ['0.org/T', 'internalField', 'uniform %g'% (TStart), ''],
        ]
    )

    # -- get the external temperature
    # cleanDataInRange = loadDataFrame(expDir, experiments[expNum])
    # bakingCurve = getTemps(cleanDataInRange, 4)
    bakingCurve = np.array([[0, 200], [10000, 200]])

    # -- change external temperature
    with open(os.path.join(baseCase.dir, "constant", "TInfTable"), "w") as fl:
        fl.writelines("(\n")
        # if not nonDeform:
        fl.writelines("\t(0\t%f)\n"%TProofing)
        fl.writelines("\t(%d\t%f)\n"%(timeProofing, TProofing))
        bakingCurve[:, 1] = TTop
        for i in range(bakingCurve.shape[0]):
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeProofing+0.1, bakingCurve[i,1]))
            fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeProofing+0.1, bakingCurve[i,1] + 273.15))
        # else:
        #     for i in range(bakingCurve.shape[0]):
        #         fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60, bakingCurve[i,1] + 273.15))
        fl.writelines(")\n")

    with open(os.path.join(baseCase.dir, "constant", "TInfTableBottom"), "w") as fl:
        fl.writelines("(\n")
        # if not nonDeform:
        fl.writelines("\t(0\t%f)\n"%TProofing)
        fl.writelines("\t(%d\t%f)\n"%(timeProofing, TProofing))
        bakingCurve[:, 1] = TBottom
        for i in range(bakingCurve.shape[0]):
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeProofing+0.1, bakingCurve[i,1]))
            # fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeProofing+0.1, 10 + bakingCurve[i,1] + 273.15))
            fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60+timeProofing+0.1, bakingCurve[i,1] + 273.15))
        # else:
        #     for i in range(bakingCurve.shape[0]):
        #         fl.write("\t(%.5g\t%.5g)\n"%(bakingCurve[i,0]*60, bakingCurve[i,1] + 273.15))
        fl.writelines(")\n")

    # 2) constant/transportProperties
    baseCase.setParameters(
        [
            ['constant/transportProperties', 'withDeformation', str(withDeformation), ''],
            ['constant/transportProperties', 'gasPermeabilityRaw', str(gasPermeabilityRaw), ''],
            ['constant/transportProperties', 'gasPermeabilityBaked', str(gasPermeabilityBaked), ''],
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
            ['system/controlDict', 'endTime', str(timeProofing), ''],
            ['system/controlDict', 'deltaT', '%.5g'%timeStepProofing, ''],
            ['system/controlDict', 'writeInterval', '%.5g'%writeIntProofing, ''],
        ]
    )

    # 6) fvSolutions
    baseCase.setParameters(
        [
            ['system/fvSolution', 'nOuterCorrectors', str(nIterProofing), 'PIMPLE'],
            ['system/fvSolution', 'D', str(DRelaxKyn), 'fields'],
            ['system/fvSolution', 'DFinal', str(DFinalRelax), 'fields'],
            ['system/fvSolution', 'omegaV', str(omegaVRelaxKyn), 'fields'],
            ['system/fvSolution', 'omegaC', str(omegaCRelaxKyn), 'fields'],
            ['system/fvSolution', 'T', str(TRelaxKyn), 'fields'],
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
            if runWithSlurm:
                baseCase.runCommands(
                    [
                        'decomposePar > log.decomposePar',
                        'srun -n%d --nodelist=%s %s -parallel > log.%s' %(nCores, node, dynSolver,dynSolver),
                    ]
                )
            else:
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
                    ['system/controlDict', 'endTime', str(timeProofing + plusTime1), ''],
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
                if runWithSlurm:
                    baseCase.runCommands(
                        [
                            'srun -n%d --nodelist=%s %s -parallel > log.%s_2' %(nCores, node, dynSolver,dynSolver),
                        ]
                    )
                else:
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
                    ['system/controlDict', 'endTime', str(timeProofing + plusTime1 + plusTime2), ''],
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
                if runWithSlurm:
                    baseCase.runCommands(
                        [
                            'srun -n%d --nodelist=%s %s -parallel > log.%s_3' %(nCores, node, dynSolver,dynSolver),
                        ]
                    )
                else:
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
            if runWithSlurm:
                baseCase.runCommands(
                    [
                        'rm -rf processor*/0',
                        'srun -n%d --nodelist=%s postProcess -parallel -func "probeOur" -dict system/probeOur > log.postProcess'%(nCores, node),
                        'srun -n%d --nodelist=%s TLFProbe -parallel -point "(1e-3 1e-3 0)" > log.TPoint6'%(nCores, node),
                        'srun -n%d --nodelist=%s intMoisture -parallel > log.intMoisture'%(nCores, node),
                        'srun -n%d --nodelist=%s getBoundPoints -parallel > log.getBoundPoints'%(nCores, node),
                    ]
                )
            else:
                baseCase.runCommands(
                    [
                        'rm -rf processor*/0',
                        'foamJob -parallel -screen postProcess -func "probeOur" -dict system/probeOur > log.postProcess',
                        'foamJob -parallel -screen TLFProbe -point "(1e-3 1e-3 0)" > log.TPoint6',
                        'foamJob -parallel -screen intMoisture > log.intMoisture',
                        'foamJob -parallel -screen getBoundPoints > log.getBoundPoints'
                    ]
                )
            for i in range(len(experiments[expNum]['probes'])):
                # thermoOffsetCorr = np.array([experiments[expNum]['thermoOffset'][0], experiments[expNum]['thermoOffset'][2], -experiments[expNum]['thermoOffset'][1]])
                thermoOffsetCorr = experiments[expNum]['thermoOffset']
                probesCorr = experiments[expNum]['probes'][i] - thermoOffsetCorr
                probesCorr[1] = np.sqrt(probesCorr[1]**2 + probesCorr[2]**2)
                probesCorr[2] = 0
                if runWithSlurm:
                    baseCase.runCommands(
                        [
                            'srun -n%d --nodelist=%s TLFProbe  -parallel -point "(%.5g %.5g %.5g)" > log.TPoint%d' %(nCores, node, probesCorr[0], probesCorr[1], probesCorr[2], i+1),
                        ]
                    )
                else:
                    baseCase.runCommands(
                        [
                            'foamJob -parallel -screen TLFProbe -point "(%.5g %.5g %.5g)" > log.TPoint%d' %(probesCorr[0], probesCorr[1], probesCorr[2], i+1),
                        ]
                    )

        if nCores > 1:
            latestTime = baseCase.latestParTime
        else:
            latestTime  = baseCase.latestTime

        shapeBef = np.loadtxt(baseCase.dir + 'Shape/D_values_%d.dat' % timeProofing, skiprows=1)   
        shapeAft = np.loadtxt(baseCase.dir + 'Shape/D_values_%d.dat' % latestTime, skiprows=1)   
        shapeBef = shapeBef[np.argsort(shapeBef[:, 0])]
        shapeAft = shapeAft[np.argsort(shapeAft[:, 0])]
        np.savetxt(os.path.join(baseCase.dir, 'Shape/D_values_%d_sorted.dat' % timeProofing), shapeBef, header='x\ty\tz', comments='')
        np.savetxt(os.path.join(baseCase.dir, 'Shape/D_values_%d_sorted.dat' % latestTime), shapeAft, header='x\ty\tz', comments='')

        saveFigPostProcess(timeProofing, outFolder)
