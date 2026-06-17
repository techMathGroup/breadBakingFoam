import os
import pandas as pd
import numpy as np

experiments = [
    # {
    #     'date': '10_02_2026',
    #     'expNumber': '1',
    #     'excelFileOven': 'opcua_log_20260210_125317.xlsx',
    #     'excelFileBread': 'Lab_test_2026_02_10.xlsm',
    # },
    {
        'date': '10_02_2026',
        'expNumber': '2',
        'excelFileOven': 'opcua_log_20260210_144124.xlsx',
        'excelFileBread': 'Lab_test_2026_02_10.xlsm',
        'stlBefore': '2026_02_10_1-Before.stl',
        'hLoaf': 4.46e-2,
        'rLoaf': 6.47e-2,
        'hLoafKynuti':4.8e-2,
        # 'rLoafKynuti': 4.4e-2,
        # 'rLoafKynuti': 4.7e-2,
        'rLoafKynuti1': 5e-2,
        'rLoafKynuti2': 5e-2,
        'up': 0,
        'probes': np.array([
            [4.45e-2, 2.89e-2, 0],
            [4.35e-2, 0, 0],
            [3.36e-2, 0, -3.63e-2],
            [1.46e-2, 0, 3.72e-2],
        ]),
        'expDispl': np.array(
            [
                [0, 0, 0],
                [6, 2.4e-2, -0.195e-2],
            ]
        ),
        # 'thermoOffset': np.array([0, 1.13, -0.31])*1e-2, 
        'thermoOffset': np.array([0, 1.13, -0.31])*1e-2, 
    },
]

def loadExpDataFromExcelMultiSheet(excelFile):
    """
    Load all sheets from Excel file into a dictionary of DataFrames.
    
    Parameters:
    -----------
    excelFile : str
        Path to Excel file
    
    Returns:
    --------
    data_dict : dict
        Dictionary with sheet names as keys and DataFrames as values
    """
    try:
        excel_file = pd.ExcelFile(excelFile)
        data_dict = {}
        
        print(f"Loading all sheets from: {excelFile}")
        for sheet_name in excel_file.sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            data_dict[sheet_name] = df
            print(f"  Sheet '{sheet_name}': {df.shape[0]} rows × {df.shape[1]} columns")
        
        return data_dict
    except Exception as e:
        print(f"Error loading Excel file: {e}")
        return None

def loadDataFrame(expDir, exp, columnWhereDataStarts=11, numberOfThermocouples=5, split=3):
    dataFromBread = loadExpDataFromExcelMultiSheet(os.path.join(expDir, exp['date'], exp['excelFileBread']))
    sheetsBread = list(dataFromBread.keys())
    expDFBread = dataFromBread[f'Data_Termo_{exp["expNumber"]}']
    # for j in range(numberOfThermocouples):
    cleanData = expDFBread[expDFBread[f'Unnamed: {columnWhereDataStarts + split*(numberOfThermocouples-1)+1}'].apply(lambda x: isinstance(x, (int, float)) and not pd.isna(x))]
    cleanDataInRange = (cleanData[cleanData[f'Unnamed: {columnWhereDataStarts + split*(numberOfThermocouples-1)+1}'] > 50])
    return cleanDataInRange

def getTemps(cleanDataInRange, j, numberOfThermocouples=5, columnWhereDataStarts=11, split=3, desiredTime=np.linspace(0, 30, 100)):        
    # for j in range(numberOfThermocouples):
    time = cleanDataInRange[f'Unnamed: {columnWhereDataStarts + j*split}'].values.astype(float) / 60
    temp = cleanDataInRange[f'Unnamed: {columnWhereDataStarts + j*split + 1}'].values.astype(float)
    tempInt = np.interp(desiredTime, time-time[0], temp)
    # print(temp, tempInt)
    return np.column_stack((desiredTime, tempInt))