import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display
import matplotlib.pyplot as plt
import seaborn as sns


def lipinski_checker(df):
    RDLogger.DisableLog('rdApp.warning')

    # 1. Create Mol objects once
    df['_mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
    valid_mask = df['_mol'].notna()

    # 2. Calculate Descriptors (using loc to avoid SettingWithCopy warnings)
    df.loc[valid_mask, 'Standard_SMILES'] = df.loc[valid_mask, '_mol'].apply(Chem.MolToSmiles)
    df.loc[valid_mask, 'MW'] = df.loc[valid_mask, '_mol'].apply(Descriptors.MolWt)
    df.loc[valid_mask, 'LogP'] = df.loc[valid_mask, '_mol'].apply(Descriptors.MolLogP)
    df.loc[valid_mask, 'HBD'] = df.loc[valid_mask, '_mol'].apply(Descriptors.NumHDonors)
    df.loc[valid_mask, 'HBA'] = df.loc[valid_mask, '_mol'].apply(Descriptors.NumHAcceptors)
    df.loc[valid_mask, 'TPSA'] = df.loc[valid_mask, '_mol'].apply(Descriptors.TPSA)
    df.loc[valid_mask, 'RotBonds'] = df.loc[valid_mask, '_mol'].apply(Descriptors.NumRotatableBonds)

    # 3. Define Lipinski Criteria (Pass = 1, Fail = 0)
    c1 = (df['MW'] <= 500).astype(int)
    c2 = (df['LogP'] <= 5).astype(int)
    c3 = (df['HBD'] <= 5).astype(int)
    c4 = (df['HBA'] <= 10).astype(int)

    # 4. Calculate Lipinski Passes (Standard Ro5 allows 1 violation)
    total_lipinski_pass = c1 + c2 + c3 + c4
    df['Lipinski_Violations'] = 4 - total_lipinski_pass
    df['Lipinski_Ro5_Passes'] = np.where(total_lipinski_pass >= 3, 'Yes', 'No')

    # 5. Define Veber Criteria
    veber_pass = (df['TPSA'] <= 140) & (df['RotBonds'] <= 10)

    # 6. DRUG-LIKENESS (Lipinski AND Veber)
    # This replaces the buggy 'if' statement
    drug_like_mask = (df['Lipinski_Ro5_Passes'] == 'Yes') & veber_pass
    df['Drug_Likeness'] = np.where(drug_like_mask, 'Yes', 'No')

    

    # Clean up
    df.drop(columns=['_mol'], inplace=True)
    RDLogger.EnableLog('rdApp.warning')
    return df

# File Handling
csv_filename = 'antibiotic_data.csv'

# Check if file exists before reading to avoid FileNotFoundError
if os.path.exists(csv_filename):
    df_data = pd.read_csv(csv_filename, index_col=0)
    RDLogger.DisableLog('rdApp.warning')
    df_data['Standard_SMILES'] = df_data['SMILES'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)) if pd.notna(x) else None)
    cleaned_df = df_data.drop_duplicates(subset=['Standard_SMILES']).copy()
    drop_duplicates = len(df_data) - len(cleaned_df)
    print(f"Number of duplicate entries found and removed: {drop_duplicates}")
    RDLogger.EnableLog('rdApp.warning')
    processed_df = lipinski_checker(cleaned_df)

    print("Processing Complete!")
    display(cleaned_df.head())
else:
    print(f"Error: File not found at {csv_filepath}")
