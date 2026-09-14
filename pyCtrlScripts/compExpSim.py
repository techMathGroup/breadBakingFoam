#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python script to post-process bread baking simulations 
# -- compare multiple experiments with multiple simulations in a several graphs
# -- legacy

#IMPORTS================================================================
import os 
import pandas as pd
import numpy as np
from OF_caseClass import *
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from expDict import *
from myAddFcs import *




# def addPlot(ax, i, j, x, y):
#     l1 = ax[i,j].plot(x, y, linewidth=2)
#     rax = plt.axes([0.05, 0.4, 0.15, 0.15])
#     labels = ['Exp 1']
#     visibility = [l1.get_visible()]
#     check = CheckButtons(rax, labels, visibility)
#     return l1

mLInit = 103.8
mSInit = 144.8

mLInit = 159.9/2
mSInit = 336.6/2

V = 7.30047E-06 * 36
alphaD0 = 0.91
alphaD0 = 0.84

nonDeform = False
# nonDeform = True

if not nonDeform:
    kynuti = 400

    simDirs = [
        '../ZZ_cases/2026/V25/exp0_nonDef_False/V88_4xDl_botZG_0.1_Close_0.1_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_18_alphaGBottom_10_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        '../ZZ_cases/2026/V25/exp1_nonDef_False/V95_2xDl_deltaConst_0.02_Close_0.02_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_15_alphaGBottom_9_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        '../ZZ_cases/2026/V25/exp2_nonDef_False/V95_2xDl_deltaConst_0.02_Close_0.02_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_15_alphaGBottom_9_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        '../ZZ_cases/2026/V25/exp3_nonDef_False/V99_2xDl_deltaConst_0.02_Close_0.02_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_7_alphaGBottom_7_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        '../ZZ_cases/2026/V25/exp4_nonDef_False/V99_2xDl_deltaConst_0.02_Close_0.02_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_7_alphaGBottom_7_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        '../ZZ_cases/2026/V25/exp5_nonDef_False/V99_2xDl_deltaConst_0.02_Close_0.02_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_70_lambda_0.42_tau_10_alphaG_7_alphaGBottom_7_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/',
        

    ]

else:
    kynuti = 0
    simDirs = [
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        '../ZZ_cases/2026/exp0_nonDef_True/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        # '../ZZ_cases/2026/exp1_nonDef_False/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        # '../ZZ_cases/2026/exp2_nonDef_False/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        # '../ZZ_cases/2026/exp3_nonDef_False/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        # '../ZZ_cases/2026/exp4_nonDef_False/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
        # '../ZZ_cases/2026/exp5_nonDef_False/mSStep_0.001_DFree_4e-06_lambda_0.55_alphaKept_0.25_tau_1_alphaG_12_alphaGBottom_20_kmsides0.008_kmbottom0.0006_kMSidesOmega0.003_kMBottomOmega_0.0006_r0_0.0009_perm_2.5e-12/',
    ]

if __name__ == "__main__":
    #INFO===============================================================
    expDir = os.path.join('..', 'Experiments2026')

    columnWhereDataStarts = 11
    numberOfThermocouples = 5
    split = 3

    fig = plt.figure(figsize=(16, 12))
    colorLst = ['r', 'g', 'b', 'm', 'c', 'y']

    #GET DATA===========================================================
    for i, exp in enumerate(experiments):

        # -- load data from oven 
        ovenExcel = os.path.join(expDir, exp['date'], f"Exp_{exp['expNumber']}", exp['excelFileOven'])
        dataFromSheets = loadExpDataFromExcelMultiSheet(ovenExcel)
        sheets = list(dataFromSheets.keys())
        expDF = dataFromSheets[sheets[0]]
        breadInOven = expDF[expDF['Weight'] > 200]
        time = breadInOven['Timestamp'].values
        weight = breadInOven['Weight'].values
        experiments[i][f'weightData'] = np.column_stack((time-time[0], weight))

        # -- load data from experimental bread
        cleanDataInRange = loadDataFrame(expDir, exp, columnWhereDataStarts, numberOfThermocouples, split)
        for j in range(numberOfThermocouples):
            experiments[i][f'data{j+1}'] = getTemps(cleanDataInRange, j, numberOfThermocouples, columnWhereDataStarts, split)

        # -- load data from simulation
        for j in range(len(exp['probes'])):
            # simDir = '../ZZ_cases/2026/exp%d/alphaG_%g_kmsides%g_kmbottom%g_kmtop%g_r0_%g_perm_%g/' % (0, 30, 4e-3, 3e-4, 3e-3, 10e-4, 1e-12)
            # simDir = '../ZZ_cases/2026/exp0_newTest/kynuti_pGFreeCoarse_newGeom_newSigma_nonMoving_nIter100_pG_tauZvednuto_tau_3_alphaG_20_kmsides0.01_kmbottom1e-10_kmtop0.003_r0_0.0001_perm_4e-13/' 
            # simDir = '../ZZ_cases/2026/exp0_newTest_5/00kyn480_pG_tS1_noOmega_tau_1_alphaG_20_kmsides0.01_kmbottom1e-10_kmtop0.003_r0_0.0003_perm_2.5e-13/' 
            simDir = simDirs[i]
            # simDir = '../ZZ_cases/2026/exp0_newTest_5/00kyn480_pG_tS1_noOmega_tau_1_alphaG_20_kmsides0.01_kmbottom1e-10_kmtop0.003_r0_0.0001_perm_2.5e-13/' 
            # simDir = '../ZZ_cases/2026/exp0_newTest_newBound/tauLow_pG_alphaG_20_kmsides0.01_kmbottom1e-10_kmtop0.003_r0_0.0004_perm_2.6e-12/' 
            TPoint = readDataFromLogFile("%s/log.TPoint%d" %(simDir, j+1))
            experiments[i][f'sim{j+1}'] = TPoint
            experiments[i][f'moisture'] = readDataFromLogFile("%s/log.intMoisture" %simDir)
            experiments[i][f'weight'] = readDataFromLogFile("%s/log.intWeigth" %simDir)

    #CREATE PLOTS=======================================================
    # Adjust spacing to make room for checkboxes
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.05, top=0.95, hspace=0.2, wspace=0.5)
    
    # Create 6 subplots
    axes = [plt.subplot(3, 2, i+1) for i in range(6)]
    
    # Store components
    all_checks = []
    
    

    for idx, ax in enumerate(axes):
        lines = []
        myLabels = []

        # TEMPERATURE===================================================
        if idx < numberOfThermocouples:
            for i in range(len(experiments)):
                # for j in range(numberOfThermocouples):
                x = experiments[i]['data'+str(idx+1)][:,0]
                y = experiments[i]['data'+str(idx+1)][:,1]
                line, = ax.plot(x, y, color=colorLst[i], linewidth=2, label= experiments[i]['date']+'_'+experiments[i]['expNumber']+ '_exp')
                myLabels.append(experiments[i]['date']+'_'+experiments[i]['expNumber']+ '_exp')
                lines.append(line)

            # Subplot styling
            ax.set_xlabel('Time (min)', fontsize=9)
            ax.set_ylabel('Temp (°C)', fontsize=9)
            ax.set_title(f'Thermocoupler {idx+1}', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
            


                
        elif idx == 5:
            for i in range(len(experiments)):
                x = experiments[i]['weightData'][:,0] 
                loss = (experiments[i]['weightData'][30,1] - experiments[i]['weightData'][:,1])
                y = (mLInit - loss) / mSInit 
                line, = ax.plot(x, y, color=colorLst[i], linewidth=2, label=experiments[i]['date']+'_'+experiments[i]['expNumber']+'_exp')
                myLabels.append(experiments[i]['date']+'_'+experiments[i]['expNumber']+'_exp')
                lines.append(line)

            for i in range(len(experiments)):
                rhoData = experiments[i][f'weight'][:,1]
                t = experiments[i][f'weight'][:,0] / 60 - kynuti / 60
                weightData = rhoData * V * alphaD0 * 1000
                loss = (weightData[0] - weightData)
                y = (mLInit - loss) / mSInit
                
                # line, = ax.plot(experiments[i][f'weight{j+1}'][:,0] / 60 - kynuti / 60, experiments[i][f'weight{j+1}'][:,1], '--',color=colorLst[i], linewidth=2, label=experiments[i]['date']+'_'+experiments[i]['expNumber']+'_sim')
                line, = ax.plot(t, y, '--',color=colorLst[i], linewidth=2, label=experiments[i]['date']+'_'+experiments[i]['expNumber']+'_sim')
                myLabels.append(experiments[i]['date']+'_'+experiments[i]['expNumber']+'_sim')
                lines.append(line)

            # Subplot styling
            ax.set_xlabel('Time (min)', fontsize=9)
            ax.set_ylabel('Weight (g)', fontsize=9)
            ax.set_title(f'Weight', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.set_ylim(0.3,0.5)
            ax.set_xlim(0, 35)

        if idx <= 3:
            for i in range(len(experiments)):
                x = experiments[i]['sim'+str(idx+1)][:,0]  / 60 - kynuti / 60
                y = experiments[i]['sim'+str(idx+1)][:,1] - 273
                line, = ax.plot(x, y, '--', color=colorLst[i], linewidth=2, label=experiments[i]['date']+'_'+experiments[i]['expNumber']+ '_sim')
                myLabels.append(experiments[i]['date']+'_'+experiments[i]['expNumber']+ '_sim')
                lines.append(line)
                ax.set_ylim(25, 110)
                ax.set_xlim(0, 35)
            
        ax.legend(ncol=2)
        
        # Calculate checkbox position
        row = idx // 2  # 0, 1, or 2
        col = idx % 2   # 0 or 1
        
        # Position checkboxes to the left of each column
        if col == 0:
            cb_x = 0  # Left column
        else:
            cb_x = 0.49  # Right column
        
        cb_y = 0.72 - row * 0.295  # Adjust for each row
        cb_width = 0.10
        cb_height = 0.2
        
        # Create checkbox axes
        rax = plt.axes([cb_x, cb_y, cb_width, cb_height])
        rax.set_title(f'Plot {idx+1}', fontsize=8)
        
        # labels = []
        # for i in range(len(experiments)):
        #     labels.append(experiments[i]['date']+'_'+experiments[i]['expNumber'])
        visibility = [True] * len(myLabels)
        check = CheckButtons(rax, myLabels, visibility)
        
        # Style checkboxes
        for rect in check.rectangles:
            rect.set_width(0.15)
            rect.set_height(0.15)
        
        # Create closure for toggle function - pass labels explicitly
        def make_toggle_func(line_list, label_list):
            def toggle(label):
                idx_map = {}
                for i, lbl in enumerate(label_list):
                    idx_map[lbl] = i
                if label in idx_map:
                    i = idx_map[label]
                    line_list[i].set_visible(not line_list[i].get_visible())
                    plt.draw()
            return toggle
        
        check.on_clicked(make_toggle_func(lines, myLabels))
        all_checks.append(check)
        

    
    plt.suptitle('Interactive Multi-Subplot Viewer - Toggle Series On/Off', 
                 fontsize=14, fontweight='bold', y=0.98)
    plt.show()
