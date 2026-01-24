#!/usr/bin/env python3
"""
Scaffold Identification and Clustering Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Load processed data
print("Loading processed data...")
df = pd.read_csv('data_with_analysis.csv')
print(f"Total compounds: {len(df)}")

# Parse molecules
df['Mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
df = df[df['Mol'].notna()].copy()

# Extract Murcko scaffolds
print("\nExtracting Murcko scaffolds...")
def get_murcko_scaffold(mol):
    """Get Murcko scaffold SMILES"""
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return None

df['Scaffold_SMILES'] = df['Mol'].apply(get_murcko_scaffold)
df = df[df['Scaffold_SMILES'].notna()].copy()

# Count scaffolds
scaffold_counts = df['Scaffold_SMILES'].value_counts()
print(f"\nTotal unique scaffolds: {len(scaffold_counts)}")
print(f"Scaffolds with >5 compounds: {len(scaffold_counts[scaffold_counts > 5])}")
print(f"Scaffolds with >10 compounds: {len(scaffold_counts[scaffold_counts > 10])}")

# Analyze scaffold distribution
print("\n" + "="*80)
print("SCAFFOLD DISTRIBUTION")
print("="*80)
print("\nTop 10 Most Common Scaffolds:")
for i, (scaffold, count) in enumerate(scaffold_counts.head(10).items(), 1):
    print(f"{i}. Scaffold: {scaffold}")
    print(f"   Count: {count} compounds")
    print(f"   Percentage: {count/len(df)*100:.1f}%")
    print()

# Assign scaffold IDs
scaffold_to_id = {scaffold: f"S{i+1}" for i, scaffold in enumerate(scaffold_counts.index)}
df['Scaffold_ID'] = df['Scaffold_SMILES'].map(scaffold_to_id)

# Save scaffold mapping
scaffold_df = pd.DataFrame({
    'Scaffold_ID': list(scaffold_to_id.values()),
    'Scaffold_SMILES': list(scaffold_to_id.keys()),
    'Compound_Count': [scaffold_counts[s] for s in scaffold_to_id.keys()]
})
scaffold_df.to_csv('scaffold_mapping.csv', index=False)
print(f"\nScaffold mapping saved to scaffold_mapping.csv")

# Analyze properties by scaffold
print("\n" + "="*80)
print("SCAFFOLD-BASED ANALYSIS")
print("="*80)

# Focus on major scaffolds (>5 compounds)
major_scaffolds = scaffold_counts[scaffold_counts > 5].index.tolist()
df_major = df[df['Scaffold_SMILES'].isin(major_scaffolds)].copy()

print(f"\nAnalyzing {len(major_scaffolds)} major scaffolds with {len(df_major)} compounds")

# Statistics by scaffold
scaffold_stats = []
for scaffold in major_scaffolds:
    scaffold_data = df[df['Scaffold_SMILES'] == scaffold].copy()
    scaffold_id = scaffold_to_id[scaffold]
    
    stats = {
        'Scaffold_ID': scaffold_id,
        'Scaffold_SMILES': scaffold,
        'N_Compounds': len(scaffold_data),
        'Mean_Sum_Ki': scaffold_data['Sum_Ki'].mean(),
        'Min_Sum_Ki': scaffold_data['Sum_Ki'].min(),
        'Max_Sum_Ki': scaffold_data['Sum_Ki'].max(),
        'Mean_cLogP': scaffold_data['cLogP'].mean(),
        'Mean_LLE': scaffold_data['Avg_LLE'].mean(),
        'Mean_LE': scaffold_data['Avg_LE'].mean(),
        'Best_Example': scaffold_data.loc[scaffold_data['Sum_Ki'].idxmin(), 'Example_Number']
    }
    scaffold_stats.append(stats)

scaffold_stats_df = pd.DataFrame(scaffold_stats)
scaffold_stats_df = scaffold_stats_df.sort_values('Mean_Sum_Ki')
scaffold_stats_df.to_csv('scaffold_statistics.csv', index=False)

print("\nScaffold Statistics (sorted by Mean Sum_Ki):")
print(scaffold_stats_df.to_string(index=False))

# Identify best scaffolds
print("\n" + "="*80)
print("BEST SCAFFOLDS FOR OPTIMIZATION")
print("="*80)

# Criteria: Low mean Sum_Ki, High mean LLE, sufficient diversity (>5 compounds)
best_scaffolds = scaffold_stats_df[
    (scaffold_stats_df['Mean_Sum_Ki'] < 0.02) & 
    (scaffold_stats_df['Mean_LLE'] > 5.0)
].copy()

print(f"\nScaffolds with Mean_Sum_Ki < 0.02 μM and Mean_LLE > 5.0: {len(best_scaffolds)}")
if len(best_scaffolds) > 0:
    print(best_scaffolds.to_string(index=False))

# Identify scaffolds worth extending
print("\n" + "="*80)
print("SCAFFOLDS WORTH EXTENDING")
print("="*80)

# Criteria: Good potency range, good LLE, sufficient compounds for SAR
extending_scaffolds = scaffold_stats_df[
    (scaffold_stats_df['Min_Sum_Ki'] < 0.01) & 
    (scaffold_stats_df['Mean_LLE'] > 4.5) &
    (scaffold_stats_df['N_Compounds'] >= 5)
].copy()

print(f"\nScaffolds with Min_Sum_Ki < 0.01 μM, Mean_LLE > 4.5, and ≥5 compounds: {len(extending_scaffolds)}")
print(extending_scaffolds.to_string(index=False))

# Identify dead-end scaffolds
print("\n" + "="*80)
print("DEAD-END SCAFFOLDS (NOT RECOMMENDED)")
print("="*80)

# Criteria: High mean Sum_Ki or low LLE
dead_end_scaffolds = scaffold_stats_df[
    (scaffold_stats_df['Mean_Sum_Ki'] > 0.1) | 
    (scaffold_stats_df['Mean_LLE'] < 3.0)
].copy()

print(f"\nScaffolds with Mean_Sum_Ki > 0.1 μM or Mean_LLE < 3.0: {len(dead_end_scaffolds)}")
if len(dead_end_scaffolds) > 0:
    print(dead_end_scaffolds[['Scaffold_ID', 'N_Compounds', 'Mean_Sum_Ki', 'Mean_cLogP', 'Mean_LLE']].to_string(index=False))

# Draw top scaffolds
print("\n" + "="*80)
print("GENERATING SCAFFOLD STRUCTURES")
print("="*80)

# Draw top 10 scaffolds
top_10_scaffolds = scaffold_stats_df.head(10)
scaffold_mols = []
scaffold_labels = []

for _, row in top_10_scaffolds.iterrows():
    mol = Chem.MolFromSmiles(row['Scaffold_SMILES'])
    if mol:
        scaffold_mols.append(mol)
        label = f"{row['Scaffold_ID']}\nN={int(row['N_Compounds'])}\nKi={row['Mean_Sum_Ki']:.4f}μM\nLLE={row['Mean_LLE']:.2f}"
        scaffold_labels.append(label)

if len(scaffold_mols) > 0:
    img = Draw.MolsToGridImage(
        scaffold_mols[:10], 
        molsPerRow=2, 
        subImgSize=(400, 400),
        legends=scaffold_labels[:10],
        returnPNG=False
    )
    img.save('top_scaffolds.png', dpi=(300, 300))
    print("Top 10 scaffolds saved to top_scaffolds.png")

# Save updated data with scaffold IDs
df.drop('Mol', axis=1).to_csv('data_with_scaffolds.csv', index=False)
print("\nData with scaffold assignments saved to data_with_scaffolds.csv")

print("\n" + "="*80)
print("Scaffold analysis complete!")
print("="*80)
