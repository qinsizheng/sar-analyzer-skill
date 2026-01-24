#!/usr/bin/env python3
"""
Generate Comprehensive Visualizations for SAR Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Load data
print("Loading data...")
df = pd.read_csv('data_with_analysis.csv')
scaffold_stats = pd.read_csv('scaffold_statistics.csv')

print(f"Total compounds: {len(df)}")

# Create output directory for figures
import os
os.makedirs('figures', exist_ok=True)

print("\n" + "="*80)
print("GENERATING VISUALIZATIONS")
print("="*80)

# 1. Potency vs Lipophilicity Plot
print("\n1. Creating Potency vs Lipophilicity plot...")
fig, ax = plt.subplots(figsize=(12, 8))

# Color by scaffold for major scaffolds
major_scaffolds = scaffold_stats[scaffold_stats['N_Compounds'] > 5]['Scaffold_ID'].tolist()
df['Scaffold_Color'] = df['Scaffold_ID'].apply(lambda x: x if x in major_scaffolds else 'Other')

# Create scatter plot
for scaffold in major_scaffolds + ['Other']:
    data = df[df['Scaffold_Color'] == scaffold]
    if len(data) > 0:
        ax.scatter(data['cLogP'], data['Sum_Ki'], 
                  label=scaffold, alpha=0.6, s=80, edgecolors='black', linewidth=0.5)

ax.set_xlabel('cLogP', fontsize=14, fontweight='bold')
ax.set_ylabel('Sum Ki (μM)', fontsize=14, fontweight='bold')
ax.set_yscale('log')
ax.set_title('BRD4 Inhibitor Potency vs Lipophilicity', fontsize=16, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)

# Add optimal region box
ax.axhline(y=0.01, color='red', linestyle='--', alpha=0.5, linewidth=2, label='Target potency')
ax.axvline(x=3.5, color='green', linestyle='--', alpha=0.5, linewidth=2, label='Target cLogP')

plt.tight_layout()
plt.savefig('figures/potency_vs_lipophilicity.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: potency_vs_lipophilicity.png")

# 2. LLE vs Sum_Ki Plot
print("\n2. Creating LLE vs Potency plot...")
fig, ax = plt.subplots(figsize=(12, 8))

for scaffold in major_scaffolds + ['Other']:
    data = df[df['Scaffold_Color'] == scaffold]
    if len(data) > 0:
        ax.scatter(data['Avg_LLE'], data['Sum_Ki'], 
                  label=scaffold, alpha=0.6, s=80, edgecolors='black', linewidth=0.5)

ax.set_xlabel('Average LLE (pKi - cLogP)', fontsize=14, fontweight='bold')
ax.set_ylabel('Sum Ki (μM)', fontsize=14, fontweight='bold')
ax.set_yscale('log')
ax.set_title('Lipophilic Ligand Efficiency vs Potency', fontsize=16, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)

# Add target regions
ax.axhline(y=0.01, color='red', linestyle='--', alpha=0.5, linewidth=2)
ax.axvline(x=6.0, color='green', linestyle='--', alpha=0.5, linewidth=2)

plt.tight_layout()
plt.savefig('figures/lle_vs_potency.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: lle_vs_potency.png")

# 3. BD1 vs BD2 Selectivity
print("\n3. Creating BD1 vs BD2 selectivity plot...")
fig, ax = plt.subplots(figsize=(10, 10))

for scaffold in major_scaffolds + ['Other']:
    data = df[df['Scaffold_Color'] == scaffold]
    if len(data) > 0:
        ax.scatter(data['Ki_BD1'], data['Ki_BD2'], 
                  label=scaffold, alpha=0.6, s=80, edgecolors='black', linewidth=0.5)

ax.set_xlabel('Ki BD1 (μM)', fontsize=14, fontweight='bold')
ax.set_ylabel('Ki BD2 (μM)', fontsize=14, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('BD1 vs BD2 Selectivity', fontsize=16, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)

# Add diagonal line (balanced)
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]
ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0, linewidth=2, label='Balanced')

plt.tight_layout()
plt.savefig('figures/bd1_vs_bd2_selectivity.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: bd1_vs_bd2_selectivity.png")

# 4. Scaffold Performance Comparison
print("\n4. Creating scaffold performance comparison...")
major_scaffold_data = scaffold_stats[scaffold_stats['N_Compounds'] > 5].copy()
major_scaffold_data = major_scaffold_data.sort_values('Mean_Sum_Ki')

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Mean Sum Ki
ax = axes[0, 0]
bars = ax.barh(major_scaffold_data['Scaffold_ID'], major_scaffold_data['Mean_Sum_Ki'], color='steelblue', edgecolor='black')
ax.set_xlabel('Mean Sum Ki (μM)', fontsize=12, fontweight='bold')
ax.set_ylabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Mean Potency by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Mean LLE
ax = axes[0, 1]
bars = ax.barh(major_scaffold_data['Scaffold_ID'], major_scaffold_data['Mean_LLE'], color='forestgreen', edgecolor='black')
ax.set_xlabel('Mean LLE', fontsize=12, fontweight='bold')
ax.set_ylabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Mean Lipophilic Ligand Efficiency by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Mean cLogP
ax = axes[1, 0]
bars = ax.barh(major_scaffold_data['Scaffold_ID'], major_scaffold_data['Mean_cLogP'], color='coral', edgecolor='black')
ax.set_xlabel('Mean cLogP', fontsize=12, fontweight='bold')
ax.set_ylabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Mean Lipophilicity by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Compound Count
ax = axes[1, 1]
bars = ax.barh(major_scaffold_data['Scaffold_ID'], major_scaffold_data['N_Compounds'], color='mediumpurple', edgecolor='black')
ax.set_xlabel('Number of Compounds', fontsize=12, fontweight='bold')
ax.set_ylabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Compound Count by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('figures/scaffold_performance_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: scaffold_performance_comparison.png")

# 5. Property Distribution Plots
print("\n5. Creating property distribution plots...")
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

properties = [
    ('MW', 'Molecular Weight (Da)', axes[0, 0]),
    ('cLogP', 'cLogP', axes[0, 1]),
    ('TPSA', 'TPSA (Ų)', axes[0, 2]),
    ('HBA', 'H-Bond Acceptors', axes[1, 0]),
    ('HBD', 'H-Bond Donors', axes[1, 1]),
    ('RotBonds', 'Rotatable Bonds', axes[1, 2])
]

for prop, label, ax in properties:
    # Create bins based on potency
    potent = df[df['Sum_Ki'] < 0.01][prop]
    moderate = df[(df['Sum_Ki'] >= 0.01) & (df['Sum_Ki'] < 0.1)][prop]
    weak = df[df['Sum_Ki'] >= 0.1][prop]
    
    ax.hist([potent, moderate, weak], bins=20, label=['Potent (<0.01 μM)', 'Moderate (0.01-0.1 μM)', 'Weak (≥0.1 μM)'],
            alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.set_xlabel(label, fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title(f'Distribution of {label}', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('figures/property_distributions.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: property_distributions.png")

# 6. Efficiency Metrics Comparison
print("\n6. Creating efficiency metrics comparison...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# LE distribution by scaffold
ax = axes[0]
scaffold_le_data = []
scaffold_labels = []
for scaffold in major_scaffolds:
    data = df[df['Scaffold_ID'] == scaffold]['Avg_LE'].dropna()
    if len(data) > 0:
        scaffold_le_data.append(data)
        scaffold_labels.append(scaffold)

bp = ax.boxplot(scaffold_le_data, labels=scaffold_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightblue')
    patch.set_edgecolor('black')
ax.set_ylabel('Average Ligand Efficiency', fontsize=12, fontweight='bold')
ax.set_xlabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Ligand Efficiency by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

# LLE distribution by scaffold
ax = axes[1]
scaffold_lle_data = []
scaffold_labels = []
for scaffold in major_scaffolds:
    data = df[df['Scaffold_ID'] == scaffold]['Avg_LLE'].dropna()
    if len(data) > 0:
        scaffold_lle_data.append(data)
        scaffold_labels.append(scaffold)

bp = ax.boxplot(scaffold_lle_data, labels=scaffold_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightgreen')
    patch.set_edgecolor('black')
ax.set_ylabel('Average LLE', fontsize=12, fontweight='bold')
ax.set_xlabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_title('Lipophilic Ligand Efficiency by Scaffold', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

plt.tight_layout()
plt.savefig('figures/efficiency_metrics_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: efficiency_metrics_comparison.png")

# 7. Top Compounds Structures
print("\n7. Drawing top compound structures...")

# Get top 12 compounds by different criteria
top_potency = df.nsmallest(12, 'Sum_Ki')
top_lle = df.nlargest(12, 'Avg_LLE')

# Draw top potency compounds
mols = []
legends = []
for idx, row in top_potency.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legend = f"Ex {int(row['Example_Number'])}\nKi={row['Sum_Ki']:.4f}μM\ncLogP={row['cLogP']:.2f}\nLLE={row['Avg_LLE']:.2f}"
        legends.append(legend)

if len(mols) > 0:
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(400, 400), 
                               legends=legends, returnPNG=False)
    img.save('figures/top_potency_compounds.png', dpi=(300, 300))
    print("   Saved: top_potency_compounds.png")

# Draw top LLE compounds
mols = []
legends = []
for idx, row in top_lle.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legend = f"Ex {int(row['Example_Number'])}\nKi={row['Sum_Ki']:.4f}μM\ncLogP={row['cLogP']:.2f}\nLLE={row['Avg_LLE']:.2f}"
        legends.append(legend)

if len(mols) > 0:
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(400, 400), 
                               legends=legends, returnPNG=False)
    img.save('figures/top_lle_compounds.png', dpi=(300, 300))
    print("   Saved: top_lle_compounds.png")

# 8. Heatmap of scaffold properties
print("\n8. Creating scaffold properties heatmap...")
major_scaffold_data = scaffold_stats[scaffold_stats['N_Compounds'] > 5].copy()
major_scaffold_data = major_scaffold_data.set_index('Scaffold_ID')

# Select key metrics for heatmap
heatmap_data = major_scaffold_data[['Mean_Sum_Ki', 'Mean_cLogP', 'Mean_LLE', 'Mean_LE', 'N_Compounds']].copy()

# Normalize for better visualization
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
heatmap_normalized = pd.DataFrame(
    scaler.fit_transform(heatmap_data),
    index=heatmap_data.index,
    columns=heatmap_data.columns
)

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(heatmap_normalized.T, annot=heatmap_data.T.round(2), fmt='g', 
            cmap='RdYlGn_r', center=0, linewidths=1, linecolor='black',
            cbar_kws={'label': 'Normalized Value'}, ax=ax)
ax.set_title('Scaffold Properties Heatmap', fontsize=16, fontweight='bold')
ax.set_xlabel('Scaffold ID', fontsize=12, fontweight='bold')
ax.set_ylabel('Property', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/scaffold_properties_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("   Saved: scaffold_properties_heatmap.png")

print("\n" + "="*80)
print("All visualizations generated successfully!")
print("="*80)
print(f"\nFigures saved in: /home/ubuntu/figures/")
