#!/usr/bin/env python3
"""
R-Group Decomposition Analysis per Scaffold
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdRGroupDecomposition as rgd
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Load data
print("Loading data with scaffold assignments...")
df = pd.read_csv('data_with_scaffolds.csv')
scaffold_stats = pd.read_csv('scaffold_statistics.csv')

# Parse molecules
df['Mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
df = df[df['Mol'].notna()].copy()

print(f"Total compounds: {len(df)}")
print(f"Total scaffolds: {df['Scaffold_ID'].nunique()}")

# Focus on major scaffolds with >5 compounds
major_scaffolds = scaffold_stats[scaffold_stats['N_Compounds'] > 5]['Scaffold_ID'].tolist()
print(f"\nAnalyzing R-groups for {len(major_scaffolds)} major scaffolds")

# R-group decomposition for each major scaffold
print("\n" + "="*80)
print("R-GROUP DECOMPOSITION ANALYSIS")
print("="*80)

all_rgroup_data = []

for scaffold_id in major_scaffolds:
    print(f"\n{'='*80}")
    print(f"Scaffold: {scaffold_id}")
    print(f"{'='*80}")
    
    # Get compounds for this scaffold
    scaffold_data = df[df['Scaffold_ID'] == scaffold_id].copy()
    scaffold_smiles = scaffold_data['Scaffold_SMILES'].iloc[0]
    
    print(f"Compounds: {len(scaffold_data)}")
    print(f"Scaffold SMILES: {scaffold_smiles}")
    
    # Get scaffold molecule
    scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
    
    if scaffold_mol is None:
        print("Could not parse scaffold molecule, skipping...")
        continue
    
    # Prepare molecules for R-group decomposition
    mols = scaffold_data['Mol'].tolist()
    
    try:
        # Perform R-group decomposition
        groups, unmatched = rgd.RGroupDecompose([scaffold_mol], mols, asSmiles=True)
        
        if len(groups) == 0:
            print("No R-groups found, trying alternative approach...")
            continue
            
        # Convert to DataFrame
        rgroup_df = pd.DataFrame(groups)
        
        # Add compound information
        rgroup_df['Example_Number'] = scaffold_data['Example_Number'].values[:len(rgroup_df)]
        rgroup_df['Ki_BD1'] = scaffold_data['Ki_BD1'].values[:len(rgroup_df)]
        rgroup_df['Ki_BD2'] = scaffold_data['Ki_BD2'].values[:len(rgroup_df)]
        rgroup_df['Sum_Ki'] = scaffold_data['Sum_Ki'].values[:len(rgroup_df)]
        rgroup_df['cLogP'] = scaffold_data['cLogP'].values[:len(rgroup_df)]
        rgroup_df['Avg_LLE'] = scaffold_data['Avg_LLE'].values[:len(rgroup_df)]
        rgroup_df['Scaffold_ID'] = scaffold_id
        
        # Identify R-group columns
        rgroup_cols = [col for col in rgroup_df.columns if col.startswith('R')]
        print(f"\nR-group positions found: {rgroup_cols}")
        
        # Analyze each R-group position
        for rgroup in rgroup_cols:
            print(f"\n{rgroup} Analysis:")
            rgroup_counts = rgroup_df[rgroup].value_counts()
            print(f"  Unique substituents: {len(rgroup_counts)}")
            
            # Show top substituents and their activity
            print(f"  Top substituents:")
            for i, (sub, count) in enumerate(rgroup_counts.head(5).items(), 1):
                sub_data = rgroup_df[rgroup_df[rgroup] == sub]
                mean_ki = sub_data['Sum_Ki'].mean()
                mean_lle = sub_data['Avg_LLE'].mean()
                print(f"    {i}. {sub} (n={count}): Mean Ki={mean_ki:.4f} μM, Mean LLE={mean_lle:.2f}")
        
        # Save R-group data
        all_rgroup_data.append(rgroup_df)
        
        # Find best and worst compounds in this scaffold
        best_idx = rgroup_df['Sum_Ki'].idxmin()
        worst_idx = rgroup_df['Sum_Ki'].idxmax()
        
        print(f"\nBest compound in {scaffold_id}:")
        print(f"  Example: {rgroup_df.loc[best_idx, 'Example_Number']}")
        print(f"  Sum_Ki: {rgroup_df.loc[best_idx, 'Sum_Ki']:.4f} μM")
        print(f"  LLE: {rgroup_df.loc[best_idx, 'Avg_LLE']:.2f}")
        for rgroup in rgroup_cols:
            print(f"  {rgroup}: {rgroup_df.loc[best_idx, rgroup]}")
        
        print(f"\nWorst compound in {scaffold_id}:")
        print(f"  Example: {rgroup_df.loc[worst_idx, 'Example_Number']}")
        print(f"  Sum_Ki: {rgroup_df.loc[worst_idx, 'Sum_Ki']:.4f} μM")
        print(f"  LLE: {rgroup_df.loc[worst_idx, 'Avg_LLE']:.2f}")
        for rgroup in rgroup_cols:
            print(f"  {rgroup}: {rgroup_df.loc[worst_idx, rgroup]}")
            
    except Exception as e:
        print(f"R-group decomposition failed: {str(e)}")
        print("Performing manual substituent analysis...")
        
        # Manual analysis: identify variable positions
        # For this dataset, we'll analyze common substitution patterns
        
        # Analyze N4 substituents (major variation point)
        print("\nManual substituent analysis:")
        
        # Sort by activity
        scaffold_data_sorted = scaffold_data.sort_values('Sum_Ki')
        
        print(f"\nTop 5 compounds:")
        for idx, row in scaffold_data_sorted.head(5).iterrows():
            print(f"  Example {row['Example_Number']}: Ki={row['Sum_Ki']:.4f} μM, cLogP={row['cLogP']:.2f}, LLE={row['Avg_LLE']:.2f}")
        
        print(f"\nBottom 5 compounds:")
        for idx, row in scaffold_data_sorted.tail(5).iterrows():
            print(f"  Example {row['Example_Number']}: Ki={row['Sum_Ki']:.4f} μM, cLogP={row['cLogP']:.2f}, LLE={row['Avg_LLE']:.2f}")

# Save combined R-group data
if len(all_rgroup_data) > 0:
    combined_rgroup_df = pd.concat(all_rgroup_data, ignore_index=True)
    combined_rgroup_df.to_csv('rgroup_decomposition.csv', index=False)
    print(f"\n\nCombined R-group data saved to rgroup_decomposition.csv")

# Perform detailed substituent analysis for the best scaffold (S6 - pyrazole series)
print("\n" + "="*80)
print("DETAILED ANALYSIS: BEST SCAFFOLD (S6 - Pyrazole Series)")
print("="*80)

s6_data = df[df['Scaffold_ID'] == 'S6'].copy()
if len(s6_data) > 0:
    print(f"\nCompounds in S6: {len(s6_data)}")
    print("\nAll S6 compounds ranked by potency:")
    s6_sorted = s6_data.sort_values('Sum_Ki')[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
    print(s6_sorted.to_string(index=False))
    
    # Analyze structural patterns
    print("\n\nKey structural features in S6:")
    # Merge back with original data to get SMILES
    s6_with_smiles = s6_data.sort_values('Sum_Ki')
    for idx, row in s6_with_smiles.iterrows():
        smiles = df.loc[idx, 'SMILES']
        # Check for key features
        has_sulfonamide = 'NS(=O)(=O)' in smiles
        has_amino = 'N' in smiles and 'NC' in smiles
        has_fluoro = 'F' in smiles
        
        features = []
        if has_sulfonamide:
            features.append('sulfonamide')
        if has_amino:
            features.append('amino')
        if has_fluoro:
            features.append('fluoro')
        
        print(f"Example {row['Example_Number']}: {', '.join(features) if features else 'base scaffold'} - Ki={row['Sum_Ki']:.4f} μM")

# Analyze the most populated scaffold (S1)
print("\n" + "="*80)
print("DETAILED ANALYSIS: MOST POPULATED SCAFFOLD (S1)")
print("="*80)

s1_data = df[df['Scaffold_ID'] == 'S1'].copy()
if len(s1_data) > 0:
    print(f"\nCompounds in S1: {len(s1_data)}")
    
    # Top 10 and bottom 10
    s1_sorted = s1_data.sort_values('Sum_Ki')
    
    print("\nTop 10 S1 compounds:")
    top10 = s1_sorted.head(10)[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
    print(top10.to_string(index=False))
    
    print("\nBottom 10 S1 compounds:")
    bottom10 = s1_sorted.tail(10)[['Example_Number', 'Chemical_Name', 'Ki_BD1', 'Ki_BD2', 'Sum_Ki', 'cLogP', 'Avg_LLE']]
    print(bottom10.to_string(index=False))
    
    # Analyze C7 substituent effects (methylsulfonylmethyl vs others)
    print("\n\nC7 Position Analysis:")
    s1_data['Has_CH2SO2CH3'] = s1_data['SMILES'].str.contains('CS\(=O\)\(=O\)C')
    
    with_ch2so2ch3 = s1_data[s1_data['Has_CH2SO2CH3'] == True]
    without_ch2so2ch3 = s1_data[s1_data['Has_CH2SO2CH3'] == False]
    
    print(f"With -CH2SO2CH3: n={len(with_ch2so2ch3)}, Mean Ki={with_ch2so2ch3['Sum_Ki'].mean():.4f} μM, Mean LLE={with_ch2so2ch3['Avg_LLE'].mean():.2f}")
    print(f"Without -CH2SO2CH3: n={len(without_ch2so2ch3)}, Mean Ki={without_ch2so2ch3['Sum_Ki'].mean():.4f} μM, Mean LLE={without_ch2so2ch3['Avg_LLE'].mean():.2f}")

print("\n" + "="*80)
print("R-group analysis complete!")
print("="*80)
