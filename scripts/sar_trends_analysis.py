#!/usr/bin/env python3
"""
SAR Trends and Efficiency Metrics Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Load data
print("Loading data...")
df = pd.read_csv('data_with_scaffolds.csv')
scaffold_stats = pd.read_csv('scaffold_statistics.csv')

print(f"Total compounds: {len(df)}")
print(f"Total scaffolds: {df['Scaffold_ID'].nunique()}")

# Focus on compounds with valid Ki data
df_valid = df[df['Sum_Ki'].notna()].copy()
print(f"Compounds with valid Ki data: {len(df_valid)}")

print("\n" + "="*80)
print("SAR TRENDS ANALYSIS")
print("="*80)

# 1. Potency vs Lipophilicity Trade-off
print("\n1. Potency vs Lipophilicity Analysis")
print("-" * 40)

# Calculate correlation
corr_ki_logp = df_valid[['Sum_Ki', 'cLogP']].corr().iloc[0, 1]
print(f"Correlation between Sum_Ki and cLogP: {corr_ki_logp:.3f}")

# Identify optimal region (high potency, low lipophilicity)
optimal_compounds = df_valid[(df_valid['Sum_Ki'] < 0.01) & (df_valid['cLogP'] < 3.5)].copy()
print(f"\nOptimal compounds (Sum_Ki < 0.01 μM, cLogP < 3.5): {len(optimal_compounds)}")
print(f"Percentage: {len(optimal_compounds)/len(df_valid)*100:.1f}%")

# Analyze by cLogP bins
df_valid['cLogP_bin'] = pd.cut(df_valid['cLogP'], bins=[-1, 1, 2, 3, 4, 10], labels=['<1', '1-2', '2-3', '3-4', '>4'])
clogp_analysis = df_valid.groupby('cLogP_bin').agg({
    'Sum_Ki': ['count', 'mean', 'median', 'min'],
    'Avg_LLE': ['mean', 'median']
}).round(4)
print("\nActivity by cLogP range:")
print(clogp_analysis)

# 2. Ligand Efficiency Analysis
print("\n\n2. Ligand Efficiency Analysis")
print("-" * 40)

# Top LE compounds
top_le = df_valid.nlargest(20, 'Avg_LE')[['Example_Number', 'Chemical_Name', 'Sum_Ki', 'cLogP', 'Avg_LE', 'Scaffold_ID']]
print("\nTop 20 compounds by Ligand Efficiency:")
print(top_le.to_string(index=False))

# LE by scaffold
le_by_scaffold = df_valid.groupby('Scaffold_ID').agg({
    'Avg_LE': ['mean', 'std', 'max'],
    'Example_Number': 'count'
}).round(4)
le_by_scaffold.columns = ['Mean_LE', 'Std_LE', 'Max_LE', 'N_Compounds']
le_by_scaffold = le_by_scaffold.sort_values('Mean_LE', ascending=False)
print("\nLigand Efficiency by Scaffold:")
print(le_by_scaffold)

# 3. Lipophilic Ligand Efficiency (LLE) Analysis
print("\n\n3. Lipophilic Ligand Efficiency (LLE) Analysis")
print("-" * 40)

# Top LLE compounds
top_lle = df_valid.nlargest(20, 'Avg_LLE')[['Example_Number', 'Chemical_Name', 'Sum_Ki', 'cLogP', 'Avg_LLE', 'Scaffold_ID']]
print("\nTop 20 compounds by Lipophilic Ligand Efficiency:")
print(top_lle.to_string(index=False))

# LLE by scaffold
lle_by_scaffold = df_valid.groupby('Scaffold_ID').agg({
    'Avg_LLE': ['mean', 'std', 'max'],
    'Example_Number': 'count'
}).round(4)
lle_by_scaffold.columns = ['Mean_LLE', 'Std_LLE', 'Max_LLE', 'N_Compounds']
lle_by_scaffold = lle_by_scaffold.sort_values('Mean_LLE', ascending=False)
print("\nLipophilic Ligand Efficiency by Scaffold:")
print(lle_by_scaffold)

# Identify compounds with LLE > 6 (excellent)
excellent_lle = df_valid[df_valid['Avg_LLE'] > 6.0].copy()
print(f"\n\nCompounds with LLE > 6.0 (excellent): {len(excellent_lle)}")
print(f"Percentage: {len(excellent_lle)/len(df_valid)*100:.1f}%")

# 4. BD1 vs BD2 Selectivity Analysis
print("\n\n4. BD1 vs BD2 Selectivity Analysis")
print("-" * 40)

# Calculate selectivity ratio
df_valid['BD1_BD2_Ratio'] = df_valid['Ki_BD1'] / df_valid['Ki_BD2']
df_valid['BD2_BD1_Ratio'] = df_valid['Ki_BD2'] / df_valid['Ki_BD1']

# Identify selective compounds (>5-fold)
bd1_selective = df_valid[df_valid['BD2_BD1_Ratio'] > 5].copy()
bd2_selective = df_valid[df_valid['BD1_BD2_Ratio'] > 5].copy()
balanced = df_valid[(df_valid['BD1_BD2_Ratio'] >= 0.2) & (df_valid['BD1_BD2_Ratio'] <= 5)].copy()

print(f"BD1-selective compounds (>5-fold): {len(bd1_selective)} ({len(bd1_selective)/len(df_valid)*100:.1f}%)")
print(f"BD2-selective compounds (>5-fold): {len(bd2_selective)} ({len(bd2_selective)/len(df_valid)*100:.1f}%)")
print(f"Balanced compounds (0.2-5 fold): {len(balanced)} ({len(balanced)/len(df_valid)*100:.1f}%)")

# Correlation between BD1 and BD2
corr_bd1_bd2 = df_valid[['Ki_BD1', 'Ki_BD2']].corr().iloc[0, 1]
print(f"\nCorrelation between Ki_BD1 and Ki_BD2: {corr_bd1_bd2:.3f}")

# Top BD2-selective compounds
if len(bd2_selective) > 0:
    top_bd2_sel = bd2_selective.nlargest(10, 'BD1_BD2_Ratio')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'BD1_BD2_Ratio', 'Scaffold_ID']]
    print("\nTop 10 BD2-selective compounds:")
    print(top_bd2_sel.to_string(index=False))

# 5. Molecular Property Analysis
print("\n\n5. Molecular Property Trends")
print("-" * 40)

# Analyze property ranges for potent compounds
potent = df_valid[df_valid['Sum_Ki'] < 0.01].copy()
moderate = df_valid[(df_valid['Sum_Ki'] >= 0.01) & (df_valid['Sum_Ki'] < 0.1)].copy()
weak = df_valid[df_valid['Sum_Ki'] >= 0.1].copy()

print(f"\nPotent compounds (Sum_Ki < 0.01 μM): {len(potent)}")
print(f"Moderate compounds (0.01 ≤ Sum_Ki < 0.1 μM): {len(moderate)}")
print(f"Weak compounds (Sum_Ki ≥ 0.1 μM): {len(weak)}")

properties = ['MW', 'cLogP', 'HBA', 'HBD', 'TPSA', 'RotBonds']
print("\nProperty comparison:")
print("\nPotent compounds:")
print(potent[properties].describe().round(2))

print("\nModerate compounds:")
print(moderate[properties].describe().round(2))

print("\nWeak compounds:")
print(weak[properties].describe().round(2))

# 6. Identify Structure-Activity Patterns
print("\n\n6. Key Structure-Activity Patterns")
print("-" * 40)

# Analyze key structural features
print("\nAnalyzing key structural motifs...")

# Methylsulfonylmethyl at C7
df_valid['Has_CH2SO2CH3'] = df_valid['SMILES'].str.contains('CS\(=O\)\(=O\)C', na=False)
with_motif = df_valid[df_valid['Has_CH2SO2CH3'] == True]
without_motif = df_valid[df_valid['Has_CH2SO2CH3'] == False]

print(f"\nWith -CH2SO2CH3 motif: n={len(with_motif)}")
print(f"  Mean Sum_Ki: {with_motif['Sum_Ki'].mean():.4f} μM")
print(f"  Median Sum_Ki: {with_motif['Sum_Ki'].median():.4f} μM")
print(f"  Mean LLE: {with_motif['Avg_LLE'].mean():.2f}")

print(f"\nWithout -CH2SO2CH3 motif: n={len(without_motif)}")
print(f"  Mean Sum_Ki: {without_motif['Sum_Ki'].mean():.4f} μM")
print(f"  Median Sum_Ki: {without_motif['Sum_Ki'].median():.4f} μM")
print(f"  Mean LLE: {without_motif['Avg_LLE'].mean():.2f}")

# Fluorine substitution
df_valid['Has_Fluorine'] = df_valid['SMILES'].str.contains('F', na=False)
with_f = df_valid[df_valid['Has_Fluorine'] == True]
without_f = df_valid[df_valid['Has_Fluorine'] == False]

print(f"\nWith fluorine: n={len(with_f)}")
print(f"  Mean Sum_Ki: {with_f['Sum_Ki'].mean():.4f} μM")
print(f"  Mean cLogP: {with_f['cLogP'].mean():.2f}")
print(f"  Mean LLE: {with_f['Avg_LLE'].mean():.2f}")

print(f"\nWithout fluorine: n={len(without_f)}")
print(f"  Mean Sum_Ki: {without_f['Sum_Ki'].mean():.4f} μM")
print(f"  Mean cLogP: {without_f['cLogP'].mean():.2f}")
print(f"  Mean LLE: {without_f['Avg_LLE'].mean():.2f}")

# 7. Best Compounds Summary
print("\n\n" + "="*80)
print("BEST COMPOUNDS SUMMARY")
print("="*80)

# Define criteria for "best" compounds
best_overall = df_valid[
    (df_valid['Sum_Ki'] < 0.005) & 
    (df_valid['cLogP'] < 3.0) & 
    (df_valid['Avg_LLE'] > 6.0)
].copy()

print(f"\nBest overall compounds (Sum_Ki < 0.005 μM, cLogP < 3.0, LLE > 6.0): {len(best_overall)}")
if len(best_overall) > 0:
    best_overall_sorted = best_overall.sort_values('Sum_Ki')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE', 'Scaffold_ID']]
    print(best_overall_sorted.to_string(index=False))

# Save analysis results
df_valid.to_csv('data_with_analysis.csv', index=False)
print("\n\nAnalysis results saved to data_with_analysis.csv")

print("\n" + "="*80)
print("SAR trends analysis complete!")
print("="*80)
