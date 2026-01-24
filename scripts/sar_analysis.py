#!/usr/bin/env python3
"""
Comprehensive SAR Analysis of BRD4 Inhibitors
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdRGroupDecomposition as rgd
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Load data
import sys
if len(sys.argv) < 2:
    print("Usage: python3 sar_analysis.py <input_csv>")
    sys.exit(1)

input_file = sys.argv[1]
print(f"Loading dataset from {input_file}...")
df = pd.read_csv(input_file)
print(f"Total compounds: {len(df)}")
print(f"Columns: {df.columns.tolist()}")

# Filter valid SMILES
df = df[df['Conversion_Success'] == True].copy()
df = df[df['SMILES'].notna()].copy()
print(f"Compounds with valid SMILES: {len(df)}")

# Parse molecules
print("\nParsing molecules...")
df['Mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
df = df[df['Mol'].notna()].copy()
print(f"Successfully parsed: {len(df)} molecules")

# Calculate molecular properties
print("\nCalculating molecular properties...")
df['MW'] = df['Mol'].apply(lambda m: Descriptors.MolWt(m))
df['cLogP'] = df['Mol'].apply(lambda m: Descriptors.MolLogP(m))
df['HBA'] = df['Mol'].apply(lambda m: Descriptors.NumHAcceptors(m))
df['HBD'] = df['Mol'].apply(lambda m: Descriptors.NumHDonors(m))
df['TPSA'] = df['Mol'].apply(lambda m: Descriptors.TPSA(m))
df['RotBonds'] = df['Mol'].apply(lambda m: Descriptors.NumRotatableBonds(m))
df['AromaticRings'] = df['Mol'].apply(lambda m: Descriptors.NumAromaticRings(m))

# Process Ki values
print("\nProcessing Ki values...")
def parse_ki(val):
    """Parse Ki values, handling > and ND cases"""
    if pd.isna(val) or val == 'ND':
        return np.nan
    val_str = str(val).strip()
    if val_str.startswith('>'):
        # For > values, use the threshold value
        return float(val_str.replace('>', '').strip())
    try:
        return float(val_str)
    except:
        return np.nan

df['Ki_BD1'] = df['Ki_BDI_uM'].apply(parse_ki)
df['Ki_BD2'] = df['Ki_BDII_uM'].apply(parse_ki)

# Calculate sum Ki and pKi values
df['Sum_Ki'] = df['Ki_BD1'] + df['Ki_BD2']
df['pKi_BD1'] = -np.log10(df['Ki_BD1'] * 1e-6)  # Convert to M
df['pKi_BD2'] = -np.log10(df['Ki_BD2'] * 1e-6)
df['Avg_pKi'] = (df['pKi_BD1'] + df['pKi_BD2']) / 2

# Calculate ligand efficiency metrics
print("\nCalculating efficiency metrics...")
df['LE_BD1'] = df['pKi_BD1'] / df['Mol'].apply(lambda m: m.GetNumHeavyAtoms())
df['LE_BD2'] = df['pKi_BD2'] / df['Mol'].apply(lambda m: m.GetNumHeavyAtoms())
df['Avg_LE'] = (df['LE_BD1'] + df['LE_BD2']) / 2

# Lipophilic Ligand Efficiency (LLE = pKi - cLogP)
df['LLE_BD1'] = df['pKi_BD1'] - df['cLogP']
df['LLE_BD2'] = df['pKi_BD2'] - df['cLogP']
df['Avg_LLE'] = (df['LLE_BD1'] + df['LLE_BD2']) / 2

# Filter compounds with valid Ki data
df_valid = df[df['Sum_Ki'].notna()].copy()
print(f"Compounds with valid Ki data: {len(df_valid)}")

# Save processed data
df_valid.drop('Mol', axis=1).to_csv('data_with_analysis.csv', index=False)
print("\nProcessed data saved to data_with_analysis.csv")

# Summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

print("\nKi Values (μM):")
print(f"Ki_BD1: {df_valid['Ki_BD1'].describe()}")
print(f"\nKi_BD2: {df_valid['Ki_BD2'].describe()}")
print(f"\nSum_Ki: {df_valid['Sum_Ki'].describe()}")

print("\nMolecular Properties:")
for prop in ['MW', 'cLogP', 'HBA', 'HBD', 'TPSA', 'RotBonds']:
    print(f"{prop}: mean={df_valid[prop].mean():.2f}, std={df_valid[prop].std():.2f}, range=[{df_valid[prop].min():.2f}, {df_valid[prop].max():.2f}]")

print("\nEfficiency Metrics:")
for metric in ['Avg_LE', 'Avg_LLE']:
    print(f"{metric}: mean={df_valid[metric].mean():.3f}, std={df_valid[metric].std():.3f}, range=[{df_valid[metric].min():.3f}, {df_valid[metric].max():.3f}]")

# Identify top compounds
print("\n" + "="*80)
print("TOP COMPOUNDS BY DIFFERENT METRICS")
print("="*80)

# Best potency (lowest Sum_Ki)
top_potency = df_valid.nsmallest(10, 'Sum_Ki')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
print("\nTop 10 Most Potent Compounds (Lowest Sum_Ki):")
print(top_potency.to_string(index=False))

# Best LLE (highest Avg_LLE)
top_lle = df_valid.nlargest(10, 'Avg_LLE')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
print("\nTop 10 Compounds by Lipophilic Ligand Efficiency (Highest Avg_LLE):")
print(top_lle.to_string(index=False))

# Best LE (highest Avg_LE)
top_le = df_valid.nlargest(10, 'Avg_LE')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LE']]
print("\nTop 10 Compounds by Ligand Efficiency (Highest Avg_LE):")
print(top_le.to_string(index=False))

# Optimal compounds: High potency + Low cLogP (Sum_Ki < 0.01 and cLogP < 4)
optimal = df_valid[(df_valid['Sum_Ki'] < 0.01) & (df_valid['cLogP'] < 4.0)].copy()
optimal = optimal.sort_values('Sum_Ki')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
print(f"\nOptimal Compounds (Sum_Ki < 0.01 μM and cLogP < 4.0): {len(optimal)} compounds")
if len(optimal) > 0:
    print(optimal.to_string(index=False))

print("\n" + "="*80)
print("Analysis complete!")
print("="*80)
