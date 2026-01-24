#!/usr/bin/env python3
"""
Generate figure with ALL unique scaffold structures from the dataset
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# Load the full dataset with scaffolds
print("Loading full dataset with scaffolds...")
df = pd.read_csv('/home/ubuntu/data_with_scaffolds.csv')

# Get unique scaffolds with their statistics
scaffold_info = df.groupby(['Scaffold_ID', 'Scaffold_SMILES']).agg({
    'Example_Number': 'count',
    'Sum_Ki': ['mean', 'min'],
    'cLogP': 'mean',
    'Avg_LLE': 'mean'
}).reset_index()

# Flatten column names
scaffold_info.columns = ['Scaffold_ID', 'Scaffold_SMILES', 'N_Compounds', 'Mean_Sum_Ki', 'Min_Sum_Ki', 'Mean_cLogP', 'Mean_LLE']

# Sort by compound count (descending)
scaffold_info = scaffold_info.sort_values('N_Compounds', ascending=False)

print(f"Total unique scaffolds: {len(scaffold_info)}")

# Prepare molecules and labels for all scaffolds
scaffold_mols = []
scaffold_labels = []

for idx, row in scaffold_info.iterrows():
    mol = Chem.MolFromSmiles(row['Scaffold_SMILES'])
    if mol:
        scaffold_mols.append(mol)
        
        # Create detailed label
        label_parts = [f"{row['Scaffold_ID']}"]
        label_parts.append(f"N={int(row['N_Compounds'])}")
        
        if pd.notna(row['Mean_Sum_Ki']):
            label_parts.append(f"Ki={row['Mean_Sum_Ki']:.3f}μM")
        if pd.notna(row['Mean_LLE']):
            label_parts.append(f"LLE={row['Mean_LLE']:.1f}")
        
        label = "\n".join(label_parts)
        scaffold_labels.append(label)

print(f"Successfully parsed {len(scaffold_mols)} scaffolds")

# Generate images in batches
print(f"\nGenerating images for {len(scaffold_mols)} scaffolds...")

batch_size = 30
num_batches = (len(scaffold_mols) + batch_size - 1) // batch_size

for batch_num in range(num_batches):
    start_idx = batch_num * batch_size
    end_idx = min((batch_num + 1) * batch_size, len(scaffold_mols))
    
    batch_mols = scaffold_mols[start_idx:end_idx]
    batch_labels = scaffold_labels[start_idx:end_idx]
    
    img = Draw.MolsToGridImage(
        batch_mols, 
        molsPerRow=5, 
        subImgSize=(350, 350),
        legends=batch_labels,
        returnPNG=False
    )
    
    filename = f'/home/ubuntu/all_scaffolds_batch{batch_num+1}.png'
    img.save(filename, dpi=(300, 300))
    print(f"Saved: all_scaffolds_batch{batch_num+1}.png (scaffolds {start_idx+1}-{end_idx})")

# Create a summary figure with top 20 most populated scaffolds
print("\nGenerating summary figure with top 20 scaffolds...")
top20_mols = scaffold_mols[:20]
top20_labels = scaffold_labels[:20]

img = Draw.MolsToGridImage(
    top20_mols, 
    molsPerRow=4, 
    subImgSize=(450, 450),
    legends=top20_labels,
    returnPNG=False
)
img.save('/home/ubuntu/top20_scaffolds_summary.png', dpi=(300, 300))
print("Saved: top20_scaffolds_summary.png")

print("\n" + "="*80)
print(f"Generated {num_batches} batch images covering all {len(scaffold_mols)} unique scaffolds!")
print("="*80)
