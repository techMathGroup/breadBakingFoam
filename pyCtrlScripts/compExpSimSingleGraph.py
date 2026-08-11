#!/usr/bin/python

# Python script to plot first experiment vs first simulation for first 4 temperatures

import os 
import pandas as pd
import numpy as np
from OF_caseClass import *
import matplotlib.pyplot as plt
from expDict import *
from myAddFcs import *

kynuti = 400
# kynuti = 1100
# kynuti = 600
# kynuti = 240
# kynuti = 300

# simDir = '../ZZ_cases/2026/V26/exp0_nonDef_False/V63_takeC_litlealphaG_Dl8e11_kMOpen_0.1_Close_0.1_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_10_lambda_0.42_tau_10_alphaG_9_alphaGBottom_9_kMSidesOmega0.01_kMBottomOmega_0.01_r0_3.3e-05_perm_1.5e-14/'
# simDir = '../ZZ_cases/2026/V26/exp0_nonDef_False/V54_testTemps__kMOpen_0.11_Close_0.11_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_10_lambda_0.42_tau_10_alphaG_9_alphaGBottom_9_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/'
# simDir = '../ZZ_cases/2026/V26/exp0_nonDef_False/V15_lowalphag_lowDl_TnonConst_0.05_Close_0.05_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_10_lambda_0.42_tau_10_alphaG_7_alphaGBottom_7_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/'
# simDir = '../ZZ_cases/2026/V19/exp0_nonDef_False/lam3KMPC_0.1_E_3000_nu_0.49_mSStep_0.0015_DFree_2e-05_tort_10_lambda_0.55_tau_10_alphaG_20_alphaGBottom_22_kMSidesOmega0.015_kMBottomOmega_0.015_r0_8e-05_perm_1e-14/'



def save_for_latex(filename, x, y, header="Time Value"):
    data = np.column_stack((x, y))
    np.savetxt(filename, data, header=header, comments='', fmt='%.6f', delimiter='\t')

def saveFigPostProcess(simDir):

    mLInit = 159.9/2
    mSInit = 336.6/2
    V = 7.30047E-06 * 36
    # alphaD0 = 0.91
    alphaD0 = 0.9

    expDir = os.path.join('..', 'Experiments2026')
    columnWhereDataStarts = 11
    numberOfThermocouples = 5
    split = 3

    # Load all experiments data
    all_experiments_data = []

    for exp in experiments:
        exp_data_dict = {}
        # Load experimental data
        cleanDataInRange = loadDataFrame(expDir, exp, columnWhereDataStarts, numberOfThermocouples, split)
        
        for j in range(4): # First 4 thermocouples
            exp_data_dict[f'data{j+1}'] = getTemps(cleanDataInRange, j, numberOfThermocouples, columnWhereDataStarts, split)

        # Load experimental moisture/weight data
        ovenExcel = os.path.join(expDir, exp['date'], f"Exp_{exp['expNumber']}", exp['excelFileOven'])
        dataFromSheets = loadExpDataFromExcelMultiSheet(ovenExcel)
        sheets = list(dataFromSheets.keys())
        expDF = dataFromSheets[sheets[0]]
        breadInOven = expDF[expDF['Weight'] > 200]
        time = breadInOven['Timestamp'].values
        weight = breadInOven['Weight'].values
        exp_data_dict['weightData'] = np.column_stack((time-time[0], weight))
        all_experiments_data.append(exp_data_dict)

    # Pick the first experiment for plotting individually
    exp = experiments[0]
    expData = all_experiments_data[0]

    # 3. Load simulation data
    simData = {}
    for j in range(4): # First 4 thermocouples
        TPoint = readDataFromLogFile("%s/log.TPoint%d" %(simDir, j+1))
        simData[f'sim{j+1}'] = TPoint
        
    simData['weight'] = readDataFromLogFile("%s/log.intWeigth" %simDir)

    # Directory to save latex data
    out_dir = 'latex_data/' + simDir.split('/')[-2]
    os.makedirs(out_dir, exist_ok=True)

    # 4. Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    colorLst = ['r', 'g', 'b', 'm']

    # --- Plot 1: Temperatures ---
    for j in range(4):
        # Experimental
        x_exp = expData[f'data{j+1}'][:,0]
        y_exp = expData[f'data{j+1}'][:,1]
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
    ax1.set_xlim(0, 35)
    ax1.legend(ncol=2)
    
    # --- Plot 2: Moisture (Weight) ---
    x_exp_w = expData['weightData'][:,0]
    loss_exp = (expData['weightData'][30,1] - expData['weightData'][:,1])
    y_exp_w = (mLInit - loss_exp) / mSInit
    ax2.plot(x_exp_w, y_exp_w, color='k', linewidth=2, label='Exp Moisture')
    save_for_latex(os.path.join(out_dir, 'moisture_exp.dat'), x_exp_w, y_exp_w, "Time(min)\tMoistureRatio")

    rhoData = simData['weight'][:,1]
    x_sim_w = simData['weight'][:,0] / 60 - kynuti / 60
    weight_sim_data = rhoData * V * alphaD0 * 1000
    print(simData['weight'])
    print(simData['weight'].shape)
    # loss_sim = (weight_sim_data[110] - weight_sim_data)
    loss_sim = (weight_sim_data[12] - weight_sim_data)
    y_sim_w = (mLInit - loss_sim) / mSInit
    ax2.plot(x_sim_w, y_sim_w, '--', color='k', linewidth=2, label='Sim Moisture')
    save_for_latex(os.path.join(out_dir, 'moisture_sim.dat'), x_sim_w, y_sim_w, "Time(min)\tMoistureRatio")

    ax2.set_xlabel('Time (min)', fontsize=12)
    ax2.set_ylabel('Moisture (dry basis)', fontsize=12)
    ax2.set_title(f"Moisture Content", fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.3, 0.5)
    ax2.set_xlim(0, 35)
    ax2.legend()
    
    # --- Create unified dataset ---
    t_unified = simData['sim1'][:,0] / 60 - kynuti / 60
    unified_data = [t_unified]
    header = ["Time(min)"]

    for j in range(4):
        x_sim = simData[f'sim{j+1}'][:,0] / 60 - kynuti / 60
        y_sim = simData[f'sim{j+1}'][:,1] - 273
        x_exp = expData[f'data{j+1}'][:,0]
        y_exp = expData[f'data{j+1}'][:,1]
        
        y_sim_interp = np.interp(t_unified, x_sim, y_sim)
        y_exp_interp = np.interp(t_unified, x_exp, y_exp)
        
        # Calculate mean and std across all experiments
        y_exp_all_interp = []
        for e in all_experiments_data:
            x_e = e[f'data{j+1}'][:,0]
            y_e = e[f'data{j+1}'][:,1]
            y_exp_all_interp.append(np.interp(t_unified, x_e, y_e))
        
        y_exp_mean = np.mean(y_exp_all_interp, axis=0)
        if j == 2:
            inds = 0,1,3,4
            # y_exp_std = np.std(y_exp_all_interp[*inds], axis=0)
            y_exp_std = np.std([y_exp_all_interp[x] for x in inds], axis=0)
        else:
            y_exp_std = np.std(y_exp_all_interp, axis=0)

        unified_data.extend([y_sim_interp, y_exp_interp, y_exp_mean, y_exp_std])
        header.extend([f"Sim_TC{j+1}", f"Exp_TC{j+1}", f"Exp_Mean_TC{j+1}", f"Exp_Std_TC{j+1}"])

    y_sim_w_interp = np.interp(t_unified, x_sim_w, y_sim_w)
    y_exp_w_interp = np.interp(t_unified, x_exp_w, y_exp_w)
    
    # Calculate mean and std across all experiments for moisture
    y_exp_w_all_interp = []
    for e in all_experiments_data:
        x_w_e = e['weightData'][:,0]
        loss_e = (e['weightData'][30,1] - e['weightData'][:,1])
        y_w_e = (mLInit - loss_e) / mSInit
        y_exp_w_all_interp.append(np.interp(t_unified, x_w_e, y_w_e))
        
    y_exp_w_mean = np.mean(y_exp_w_all_interp, axis=0)
    # print(y_exp_w_all_interp)
    # print(y_exp_w_all_interp.shape())
    y_exp_w_std = np.std(y_exp_w_all_interp[0:4], axis=0)
    
    unified_data.extend([y_sim_w_interp, y_exp_w_interp, y_exp_w_mean, y_exp_w_std])
    header.extend(["Sim_Moisture", "Exp_Moisture", "Exp_Mean_Moisture", "Exp_Std_Moisture"])
    
    unified_data_array = np.column_stack(unified_data)
    np.savetxt(os.path.join(out_dir, 'unified_data.dat'), unified_data_array, header='\t'.join(header), comments='', fmt='%.6f', delimiter='\t')

    plt.suptitle(f"Experiment {exp['date']}_{exp['expNumber']} vs Simulation", fontsize=16, fontweight='bold')
    plt.subplots_adjust(top=0.9)

    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join(simDir, 'exp_vs_sim.png'))


# simDir = '../ZZ_cases/2026/V26/exp0_nonDef_False/V14_lowalphag_lowDl_TConst_0.05_Close_0.05_E_3000_nu_0.49_mSStep_0.001_DFree_2.6e-05_tortOpen_3_tortClosed_10_lambda_0.42_tau_10_alphaG_7_alphaGBottom_7_kMSidesOmega0.01_kMBottomOmega_0.01_r0_1.1e-05_perm_1.5e-15/'
# saveFigPostProcess(simDir)