#!/usr/bin/env python3
"""
Create SAR Tables for S1 (most populated) and S7 (with all R-groups)
"""

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import warnings
warnings.filterwarnings('ignore')

# Load data
print("Loading R-group decomposition data...")
df_all = pd.read_csv('rgroup_decomposition.csv')

# Process S1 and S7
scaffolds = ['S1', 'S7']

for scaffold_id in scaffolds:
    print(f"\n{'='*80}")
    print(f"Processing Scaffold {scaffold_id}")
    print(f"{'='*80}")
    
    # Filter data for this scaffold
    scaffold_data = df_all[df_all['Scaffold_ID'] == scaffold_id].copy()
    
    if len(scaffold_data) == 0:
        print(f"No data found for {scaffold_id}, skipping...")
        continue
    
    scaffold_data = scaffold_data.sort_values('Sum_Ki')
    print(f"Compounds: {len(scaffold_data)}")
    
    # Get core structure
    core_smiles = scaffold_data['Core'].iloc[0]
    print(f"Core SMILES: {core_smiles}")
    
    # For S7, include ALL R-groups (R1-R9)
    if scaffold_id == 'S7':
        rgroup_cols = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']
    else:
        # For S1, get all R-group columns
        rgroup_cols = [col for col in scaffold_data.columns if col.startswith('R')]
    
    # Filter out completely empty R-groups
    rgroup_cols = [col for col in rgroup_cols if col in scaffold_data.columns and scaffold_data[col].notna().any()]
    print(f"R-group positions: {rgroup_cols}")
    
    # Generate core structure image
    core_mol = Chem.MolFromSmiles(core_smiles)
    if core_mol:
        img_core = Draw.MolToImage(core_mol, size=(1400, 1200), kekulize=True)
        img_core.save(f'{scaffold_id}_Core_Structure.png', dpi=(600, 600))
        print(f"Saved: {scaffold_id}_Core_Structure.png")
    
    # Prepare table data
    table_data = []
    for idx, row in scaffold_data.iterrows():
        row_data = [int(row['Example_Number'])]
        
        # Add R-groups
        for rgroup in rgroup_cols:
            if rgroup in row and pd.notna(row[rgroup]):
                r_val = str(row[rgroup]).replace('[H]', 'H')
                # Clean up the R-group notation
                for i in range(1, 10):
                    r_val = r_val.replace(f'[*:{i}]', '')
                row_data.append(r_val)
            else:
                row_data.append('-')
        
        # Add activity data
        ki_bd1_nm = row['Ki_BD1'] * 1000
        ki_bd2_nm = row['Ki_BD2'] * 1000
        sum_ki_nm = row['Sum_Ki'] * 1000
        clogp = row['cLogP']
        lle = row['Avg_LLE']
        
        row_data.extend([
            f"{ki_bd1_nm:.2f}",
            f"{ki_bd2_nm:.2f}",
            f"{sum_ki_nm:.2f}",
            f"{clogp:.2f}",
            f"{lle:.2f}"
        ])
        
        table_data.append(row_data)
    
    # Create figure - larger for S1 due to many compounds
    if scaffold_id == 'S1':
        fig_height = max(20, 8 + len(scaffold_data) * 0.5)
    else:
        fig_height = max(16, 8 + len(scaffold_data) * 0.8)
    
    fig = plt.figure(figsize=(28, fig_height))
    
    # Title
    title_text = f'SAR Table for Scaffold {scaffold_id} - BRD4 Inhibitors'
    if scaffold_id == 'S1':
        title_text += ' (Most Populated Series)'
    fig.text(0.5, 0.98, title_text, ha='center', va='top', fontsize=24, fontweight='bold')
    
    # Add core structure
    if core_mol:
        ax_core = fig.add_axes([0.35, 0.90, 0.3, 0.06])
        ax_core.imshow(img_core)
        ax_core.axis('off')
    
    # Create column headers
    columns = ['Cmpd\n#']
    for rgroup in rgroup_cols:
        columns.append(f'{rgroup}')
    columns.extend(['Ki BD1\n(nM)', 'Ki BD2\n(nM)', 'Sum Ki\n(nM)', 'cLogP', 'LLE'])
    
    # Calculate table position based on number of rows
    if scaffold_id == 'S1':
        table_height = 0.80
        table_y = 0.02
    else:
        table_height = min(0.75, 0.15 + len(scaffold_data) * 0.05)
        table_y = 0.05
    
    ax_table = fig.add_axes([0.02, table_y, 0.96, table_height])
    ax_table.axis('tight')
    ax_table.axis('off')
    
    # Calculate column widths dynamically
    n_rgroups = len(rgroup_cols)
    n_total_cols = n_rgroups + 6  # R-groups + Cmpd# + 5 activity columns
    
    # Allocate widths
    cmpd_width = 0.04
    activity_width = 0.06
    rgroup_total = 1.0 - cmpd_width - (5 * activity_width)
    rgroup_width = rgroup_total / n_rgroups
    
    col_widths = [cmpd_width] + [rgroup_width] * n_rgroups + [activity_width] * 5
    
    # Create table
    table = ax_table.table(cellText=table_data, colLabels=columns,
                          cellLoc='center', loc='center',
                          colWidths=col_widths)
    
    table.auto_set_font_size(False)
    if scaffold_id == 'S1':
        table.set_fontsize(8)
        table.scale(1, 2.0)
    else:
        table.set_fontsize(10)
        table.scale(1, 3.0)
    
    # Style header
    for i in range(len(columns)):
        cell = table[(0, i)]
        cell.set_facecolor('#2E5090')
        cell.set_text_props(weight='bold', color='white')
        cell.set_edgecolor('white')
        cell.set_linewidth(2)
    
    # Style data rows
    for i in range(1, len(table_data) + 1):
        for j in range(len(columns)):
            cell = table[(i, j)]
            cell.set_edgecolor('gray')
            cell.set_linewidth(0.5)
            
            if i % 2 == 0:
                cell.set_facecolor('#F2F2F2')
            else:
                cell.set_facecolor('white')
            
            # Highlight top 3 compounds
            if i <= 3:
                if i == 1:
                    cell.set_facecolor('#FFD700')  # Gold for best
                elif i == 2:
                    cell.set_facecolor('#FFA500')  # Orange for 2nd
                elif i == 3:
                    cell.set_facecolor('#FFE4B5')  # Light orange for 3rd
                cell.set_text_props(weight='bold')
            
            # Highlight excellent Sum Ki values
            if j == len(columns) - 3:  # Sum Ki column
                try:
                    val = float(table_data[i-1][j])
                    if val < 5 and i > 3:  # Don't re-highlight top 3
                        cell.set_facecolor('#90EE90')
                        cell.set_text_props(weight='bold', color='darkgreen')
                except:
                    pass
            
            # Highlight excellent LLE values
            if j == len(columns) - 1:  # LLE column
                try:
                    val = float(table_data[i-1][j])
                    if val > 7.0 and i > 3:  # Don't re-highlight top 3
                        cell.set_facecolor('#90EE90')
                        cell.set_text_props(weight='bold', color='darkgreen')
                except:
                    pass
    
    # Add summary statistics
    mean_ki = scaffold_data['Sum_Ki'].mean() * 1000
    best_ki = scaffold_data['Sum_Ki'].min() * 1000
    mean_lle = scaffold_data['Avg_LLE'].mean()
    top3_ki = scaffold_data.nsmallest(3, 'Sum_Ki')['Sum_Ki'].tolist()
    
    summary_text = f"""
    Summary Statistics for {scaffold_id}:
    • N compounds: {len(scaffold_data)}
    • Best Sum Ki: {best_ki:.2f} nM (Example {int(scaffold_data.iloc[0]['Example_Number'])})
    • Top 3 Sum Ki: {', '.join([f'{k*1000:.2f}' for k in top3_ki])} nM
    • Mean Sum Ki: {mean_ki:.2f} nM
    • Mean LLE: {mean_lle:.2f}
    
    Color coding:
    • Gold = Best | Orange = 2nd | Light Orange = 3rd
    • Green = Excellent values (Sum Ki < 5 nM or LLE > 7.0)
    """
    
    fig.text(0.98, 0.88, summary_text, fontsize=11, verticalalignment='top', ha='right',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, pad=1))
    
    # Save figure
    plt.savefig(f'SAR_Table_{scaffold_id}_Complete.png', 
               dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: SAR_Table_{scaffold_id}_Complete.png")
    plt.close()

print("\n" + "="*80)
print("S1 AND S7 SAR TABLES GENERATED SUCCESSFULLY!")
print("="*80)
print("\nGenerated files:")
print("  - SAR_Table_S1_Complete.png (55 compounds)")
print("  - S1_Core_Structure.png")
print("  - SAR_Table_S7_Complete.png (all R1-R9 groups)")
print("  - S7_Core_Structure.png")
