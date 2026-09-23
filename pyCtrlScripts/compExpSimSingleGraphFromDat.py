#!/usr/bin/python

# Python script to plot compare experimental temperatures with simulation

import os 
import pandas as pd
import numpy as np
from OF_caseClass import *
import matplotlib.pyplot as plt
from expDict import *
from myAddFcs import *


def save_for_latex(filename, x, y, header="Time Value"):
    data = np.column_stack((x, y))
    np.savetxt(filename, data, header=header, comments='', fmt='%.6f', delimiter='\t')

def saveFigPostProcess(kynuti, simDir):
    expProbeInds = [3, 7, 11, 15] 
    expMois = 18
    mLInit = 159.9/2
    mSInit = 336.6/2
    # V = 7.30047E-06 * 36
    # alphaD0 = 0.91
    # alphaD0 = 0.9
    alphaD0 = 0.84

    expDir = os.path.join('..', 'Experiments2026')
    columnWhereDataStarts = 11
    numberOfThermocouples = 5
    split = 3

    # Load all experiments data
    all_experiments_data = []

    for exp in experiments:
        exp_data_dict = {}
        # Load experimental data
        # cleanDataInRange = loadDataFrame(expDir, exp, columnWhereDataStarts, numberOfThermocouples, split)

        expData = np.loadtxt(os.path.join(simDir, "ZZ_dataForPostProcessing", "V2.dat"), skiprows=1)
        
        # for j in range(4): # First 4 thermocouples
        #     exp_data_dict[f'data{j+1}'] = 


        # Load experimental moisture/weight data
        # ovenExcel = os.path.join(expDir, exp['date'], f"Exp_{exp['expNumber']}", exp['excelFileOven'])
        # dataFromSheets = loadExpDataFromExcelMultiSheet(ovenExcel)
        # sheets = list(dataFromSheets.keys())
        # expDF = dataFromSheets[sheets[0]]
        # breadInOven = expDF[expDF['Weight'] > 200]
        # time = breadInOven['Timestamp'].values
        # weight = breadInOven['Weight'].values
        # exp_data_dict['weightData'] = np.column_stack((time-time[0], weight))
        # all_experiments_data.append(exp_data_dict)

    # Pick the first experiment for plotting individually
    exp = experiments[0]
    # expData = all_experiments_data[0]

    # 3. Load simulation data
    simData = {}
    for j in range(4): # First 4 thermocouples
        TPoint = readDataFromLogFile("%s/log.TPoint%d" %(simDir, j+1))
        simData[f'sim{j+1}'] = TPoint
        
    # simData['weight'] = readDataFromLogFile("%s/log.intWeigth" %simDir)
    simData['moisture'] = readDataFromLogFile("%s/log.intMoisture" %simDir)


    # Directory to save latex data
    out_dir = 'latex_data/' + simDir.split('/')[-2]
    os.makedirs(out_dir, exist_ok=True)

    # 4. Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    colorLst = ['r', 'g', 'b', 'm']

    # --- Plot 1: Temperatures ---
    for j in range(4):
        # Experimental
        x_exp = expData[:,0]
        y_exp = expData[:, expProbeInds[j]]
        ax1.plot(x_exp, y_exp, color=colorLst[j], linewidth=2, label=f'Exp TC{j+1}')
        save_for_latex(os.path.join(out_dir, f'temp_exp_TC{j+1}.dat'), x_exp, y_exp, "Time(min)\tTemp(C)")
        
        # Simulation
        x_sim = simData[f'sim{j+1}'][:,0] / 60 - kynuti / 60
        y_sim = simData[f'sim{j+1}'][:,1] - 273
        ax1.plot(x_sim, y_sim, '--', color=colorLst[j], linewidth=2, label=f'Sim TC{j+1}')
        save_for_latex(os.path.join(out_dir, f'temp_sim_TC{j+1}.dat'), x_sim, y_sim, "Time(min)\tTemp(C)")

    ax1.set_xlabel('Time (min)', fontsize=12)
    ax1.set_ylabel('Temp (°C)', fontsize=12)
    ax1.set_title(f"Temperatures", fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(25, 110)
    ax1.set_xlim(0, 25)
    ax1.legend(ncol=2)
    
    # --- Plot 2: Moisture (Weight) ---
    x_exp_w = expData[:,0]
    # loss_exp = (expData['weightData'][30,1] - expData['weightData'][:,1])
    y_exp_w = expData[:, expMois]
    ax2.plot(x_exp_w, y_exp_w, color='k', linewidth=2, label='Exp Moisture')
    save_for_latex(os.path.join(out_dir, 'moisture_exp.dat'), x_exp_w, y_exp_w, "Time(min)\tMoistureRatio")

    # rhoData = simData['weight'][:,1]
    # x_sim_w = simData['weight'][:,0] / 60 - kynuti / 60
    # weight_sim_data = rhoData * 36  * 1000
    # print(simData['weight'])
    # print(simData['weight'].shape)
    # print(simData['moisture'])
    # print(simData['moisture'].shape)
    # loss_sim = (weight_sim_data[110] - weight_sim_data)
    # loss_sim = (weight_sim_data[10] - weight_sim_data)
    # loss_sim = (weight_sim_data[0] - weight_sim_data)
    # y_sim_w = rhoData * 72 * 1000
    # ax3.plot(x_sim_w, y_sim_w, '--', color='k', linewidth=2, label='Sim Moisture')
    # ax3.plot(expData['weightData'][:,0], expData['weightData'][:,1]-expData['weightData'][5,1]+y_sim_w[0], '--', color='k', linewidth=2, label='Sim Moisture')
    ax2.plot(simData['moisture'][:,0] / 60 - kynuti / 60, simData['moisture'][:,1], '--', color='k', linewidth=2, label='Sim Moisture')
    # save_for_latex(os.path.join(out_dir, 'moisture_sim.dat'), x_sim_w, y_sim_w, "Time(min)\tMoistureRatio")
    # ax3.set_xlim(0, 25)

    ax2.set_xlabel('Time (min)', fontsize=12)
    ax2.set_ylabel('Moisture (dry basis)', fontsize=12)
    ax2.set_title(f"Moisture Content", fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.3, 0.5)
    ax2.set_xlim(0, 25)
    ax2.legend()
    
    # --- Create unified dataset ---
    t_unified = simData['sim1'][:,0] / 60 - kynuti / 60
    unified_data = [t_unified]
    header = ["Time(min)"]

    # x_oven = expData[f'data0'][:,0]
    # y_oven = expData[f'data0'][:,1]
    # y_oven = np.interp(t_unified, x_oven, y_oven)

    plt.suptitle(f"Experiment {exp['date']}_{exp['expNumber']} vs Simulation", fontsize=16, fontweight='bold')
    plt.subplots_adjust(top=0.9)

    plt.tight_layout()
    # plt.show()
    # plt.savefig(os.path.join(simDir, 'exp_vs_sim3.png'))
    plt.savefig(os.path.join(simDir, 'exp_vs_sim4.png'))
    # plt.savefig(os.path.join(simDir, 'exp_vs_sim2.png'))

# saveFigPostProcess(2400, "../ZZ_cases/01_breadAx2DOurExp2/")