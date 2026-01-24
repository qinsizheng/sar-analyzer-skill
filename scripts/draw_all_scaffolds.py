#!/usr/bin/env python3
"""
Generate figure with all scaffold structures
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# Load scaffold data
print("Loading scaffold data...")
scaffold_stats = pd.read_csv('/home/ubuntu/scaffold_statistics.csv')

# Sort by compound count (descending) to show most important scaffolds first
scaffold_stats = scaffold_stats.sort_values('N_Compounds', ascending=False)

print(f"Total scaffolds: {len(scaffold_stats)}")

# Prepare molecules and labels for all scaffolds
scaffold_mols = []
scaffold_labels = []

for idx, row in scaffold_stats.iterrows():
    mol = Chem.MolFromSmiles(row['Scaffold_SMILES'])
    if mol:
        scaffold_mols.append(mol)
        
        # Create detailed label
        label_parts = [f"{row['Scaffold_ID']}"]
        label_parts.append(f"N={int(row['N_Compounds'])}")
        
        # Add statistics
        label_parts.append(f"Ki={row['Mean_Sum_Ki']:.4f}μM")
        label_parts.append(f"LLE={row['Mean_LLE']:.2f}")
        
        label = "\n".join(label_parts)
        scaffold_labels.append(label)

print(f"Successfully parsed {len(scaffold_mols)} scaffolds")

# Generate comprehensive grid image
# For large numbers of scaffolds, we'll create multiple images or use a large grid
if len(scaffold_mols) <= 30:
    # Single image with all scaffolds
    print("\nGenerating single image with all scaffolds...")
    img = Draw.MolsToGridImage(
        scaffold_mols, 
        molsPerRow=5, 
        subImgSize=(400, 400),
        legends=scaffold_labels,
        returnPNG=False
    )
    img.save('/home/ubuntu/all_scaffolds_complete.png', dpi=(300, 300))
    print("Saved: all_scaffolds_complete.png")
    
else:
    # Split into multiple images
    print(f"\nGenerating multiple images for {len(scaffold_mols)} scaffolds...")
    
    # First image: Top 30 most populated scaffolds
    img = Draw.MolsToGridImage(
        scaffold_mols[:30], 
        molsPerRow=5, 
        subImgSize=(400, 400),
        legends=scaffold_labels[:30],
        returnPNG=False
    )
    img.save('/home/ubuntu/all_scaffolds_top30.png', dpi=(300, 300))
    print("Saved: all_scaffolds_top30.png (top 30 by compound count)")
    
    # Second image: Next 30 scaffolds
    if len(scaffold_mols) > 30:
        batch2_mols = scaffold_mols[30:60]
        batch2_labels = scaffold_labels[30:60]
        
        img = Draw.MolsToGridImage(
            batch2_mols, 
            molsPerRow=5, 
            subImgSize=(400, 400),
            legends=batch2_labels,
            returnPNG=False
        )
        img.save('/home/ubuntu/all_scaffolds_batch2.png', dpi=(300, 300))
        print(f"Saved: all_scaffolds_batch2.png (scaffolds 31-60)")
    
    # Third image: Next 30 scaffolds
    if len(scaffold_mols) > 60:
        batch3_mols = scaffold_mols[60:90]
        batch3_labels = scaffold_labels[60:90]
        
        img = Draw.MolsToGridImage(
            batch3_mols, 
            molsPerRow=5, 
            subImgSize=(400, 400),
            legends=batch3_labels,
            returnPNG=False
        )
        img.save('/home/ubuntu/all_scaffolds_batch3.png', dpi=(300, 300))
        print(f"Saved: all_scaffolds_batch3.png (scaffolds 61-90)")
    
    # Fourth image: Remaining scaffolds
    if len(scaffold_mols) > 90:
        batch4_mols = scaffold_mols[90:]
        batch4_labels = scaffold_labels[90:]
        
        img = Draw.MolsToGridImage(
            batch4_mols, 
            molsPerRow=5, 
            subImgSize=(400, 400),
            legends=batch4_labels,
            returnPNG=False
        )
        img.save('/home/ubuntu/all_scaffolds_batch4.png', dpi=(300, 300))
        print(f"Saved: all_scaffolds_batch4.png (remaining {len(batch4_mols)} scaffolds)")

# Also create a focused image of just the major scaffolds (>5 compounds)
major_scaffolds = scaffold_stats[scaffold_stats['N_Compounds'] > 5]
print(f"\nGenerating focused image for {len(major_scaffolds)} major scaffolds...")

major_mols = []
major_labels = []

for idx, row in major_scaffolds.iterrows():
    mol = Chem.MolFromSmiles(row['Scaffold_SMILES'])
    if mol:
        major_mols.append(mol)
        
        label_parts = [f"{row['Scaffold_ID']}"]
        label_parts.append(f"N={int(row['N_Compounds'])}")
        label_parts.append(f"Mean Ki={row['Mean_Sum_Ki']:.4f}μM")
        label_parts.append(f"Best Ki={row['Min_Sum_Ki']:.4f}μM")
        label_parts.append(f"LLE={row['Mean_LLE']:.2f}")
        
        label = "\n".join(label_parts)
        major_labels.append(label)

if len(major_mols) > 0:
    img = Draw.MolsToGridImage(
        major_mols, 
        molsPerRow=3, 
        subImgSize=(500, 500),
        legends=major_labels,
        returnPNG=False
    )
    img.save('/home/ubuntu/major_scaffolds_detailed.png', dpi=(300, 300))
    print("Saved: major_scaffolds_detailed.png")

print("\n" + "="*80)
print("All scaffold figures generated successfully!")
print("="*80)
