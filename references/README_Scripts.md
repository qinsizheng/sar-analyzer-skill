# SAR Analysis Scripts Package

## Overview
This package contains all Python scripts and documentation used for the comprehensive Structure-Activity Relationship (SAR) analysis of 276 BRD4 inhibitors.

## Documentation Files

1. **SAR_Analysis_Methods.md**
   - Complete methodology documentation
   - Detailed description of all analysis steps
   - References to specific scripts
   - Key findings summary

2. **LLE_Calculation_Methodology.md**
   - Detailed explanation of LLE calculation
   - Step-by-step examples
   - Interpretation guidelines

## Core Analysis Scripts

### Phase 1: Data Loading and Property Calculation
**sar_analysis.py**
- Load CSV data with SMILES structures
- Parse Ki values (handle <, >, ≥ symbols)
- Calculate molecular properties (MW, cLogP, HBA, HBD, TPSA)
- Compute efficiency metrics (LE, LLE)
- Output: data_with_analysis.csv

### Phase 2: Scaffold Identification
**scaffold_analysis.py**
- Generate Murcko scaffolds
- Cluster compounds by scaffold
- Identify major scaffolds (≥5 compounds)
- Calculate scaffold statistics
- Outputs: scaffold_statistics.csv, scaffold_mapping.csv, data_with_scaffolds.csv

### Phase 3: R-Group Decomposition
**rgroup_analysis.py**
- Define core structures with attachment points
- Perform R-group decomposition using RDKit
- Extract and clean R-group SMILES
- Map R-groups to activity data
- Output: rgroup_decomposition.csv

### Phase 4: Statistical Analysis
**sar_trends_analysis.py**
- Correlation analysis (Pearson, Spearman)
- Scaffold comparison (ANOVA)
- Identify top compounds by different metrics
- Generate summary statistics
- Output: Console statistics

### Phase 5: Visualization
**create_visualizations.py**
- Potency vs lipophilicity scatter plots
- Scaffold comparison bar charts
- Top compounds visualization
- Export high-resolution PNG files
- Outputs: Multiple figures in /figures/ directory

## SAR Table Generation Scripts

### Scaffold-Specific Tables
**create_s1_s7_tables.py**
- Generate detailed SAR tables for S1 (55 compounds) and S7 (all R1-R9)
- Core structure with R-group labels
- Activity data with color-coded highlighting
- Summary statistics
- Outputs: SAR_Table_S1_Complete.png, SAR_Table_S7_Complete.png

### Scaffold Structure Visualization
**draw_all_scaffolds.py**
- Generate images of major scaffolds (S1-S7)
- High-resolution structure images
- Statistics labels
- Outputs: S*_Core_Structure.png files

**draw_all_unique_scaffolds.py**
- Generate comprehensive figures for all 113 unique scaffolds
- Batch images (30 scaffolds per image)
- Top 20 summary figure
- Outputs: all_scaffolds_batch*.png, top20_scaffolds_summary.png

## Dependencies

Required Python packages:
- rdkit (2023.9.1 or later)
- pandas (2.1.0 or later)
- numpy (1.24.3 or later)
- matplotlib (3.7.2 or later)
- seaborn (0.12.2 or later)
- scipy (1.11.2 or later)

Install with:
```bash
pip3 install rdkit pandas numpy matplotlib seaborn scipy
```

## Usage Instructions

### Step 1: Prepare Input Data
Ensure you have the input CSV file with columns:
- Example_Number
- Chemical_Name
- SMILES
- Ki_BDI_uM
- Ki_BDII_uM

### Step 2: Run Core Analysis (in order)
```bash
python3 sar_analysis.py
python3 scaffold_analysis.py
python3 rgroup_analysis.py
python3 sar_trends_analysis.py
python3 create_visualizations.py
```

### Step 3: Generate SAR Tables
```bash
python3 create_s1_s7_tables.py
python3 draw_all_scaffolds.py
python3 draw_all_unique_scaffolds.py
```

## Output Files

### Data Files
- data_with_analysis.csv: Full dataset with calculated properties
- scaffold_statistics.csv: Scaffold-level summary
- scaffold_mapping.csv: Compound-to-scaffold mapping
- rgroup_decomposition.csv: R-group analysis results

### Figure Files
- potency_vs_lipophilicity.png
- scaffold_comparison.png
- top_potency_compounds.png
- SAR_Table_S1_Complete.png
- SAR_Table_S7_Complete.png
- all_scaffolds_batch*.png
- top20_scaffolds_summary.png

## Key Findings

### Best Scaffolds
1. S6 (Pyrazole): Best LLE (6.35), excellent potency
2. S3 (Cyclohexyl): High LLE (6.37), consistent activity
3. S1 (Parent): Most explored, proven track record

### Top Compounds
1. Example 9 (S1): Sum Ki = 1.56 nM, LLE = 7.20
2. Example 22 (S6): Sum Ki = 1.82 nM, LLE = 7.27
3. Example 24 (S2): Sum Ki = 1.76 nM, LLE = 7.16

### SAR Insights
- R1: CS(=O)(=O)C or CCS(=O)(=O)NH optimal
- R2: H preferred; bulky groups reduce activity
- R3: F (para) on phenyl ring beneficial
- Keep cLogP < 3.0 for optimal LLE

## Version Information

- Analysis Date: January 15, 2026
- Python Version: 3.11.0
- RDKit Version: 2023.9.1
- Total Compounds: 276
- Major Scaffolds: 7 (S1-S7)
