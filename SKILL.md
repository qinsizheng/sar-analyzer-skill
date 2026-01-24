---
name: sar-analyzer
description: >
  Comprehensive Structure-Activity Relationship (SAR) analysis for small molecule drug discovery.
  Performs scaffold identification, R-group decomposition, molecular property calculation,
  and efficiency metrics (LE, LLE) for medicinal chemistry optimization.
  Supports analysis of compound libraries with Ki/IC50 data, generates publication-quality
  SAR tables with structures, and identifies optimal scaffolds for lead optimization.
  Keywords: cheminformatics, drug-discovery, SAR, medicinal-chemistry, RDKit, scaffold-analysis, R-group-decomposition
---

# SAR Analyzer Skill

## Overview

The **SAR Analyzer** skill provides a complete computational workflow for Structure-Activity Relationship (SAR) analysis of small molecule libraries. It is designed for medicinal chemists and computational chemists working on drug discovery projects who need to:

- **Identify and cluster molecular scaffolds** using Murcko framework extraction
- **Perform R-group decomposition** to understand substituent contributions to activity
- **Calculate molecular properties** and efficiency metrics (LE, LLE)
- **Generate publication-quality SAR tables** with chemical structures
- **Visualize potency-lipophilicity relationships** and scaffold comparisons
- **Identify optimal compounds** and prioritize scaffolds for further optimization

This skill integrates RDKit-based cheminformatics with statistical analysis and visualization to provide actionable insights for medicinal chemistry campaigns.

---

## When to Use This Skill

Use the **SAR Analyzer** skill when the user requests:

- SAR analysis of a compound library
- Scaffold identification or clustering
- R-group decomposition or substituent analysis
- Molecular property calculations (MW, cLogP, TPSA, etc.)
- Ligand efficiency (LE) or lipophilic ligand efficiency (LLE) calculations
- Generation of SAR tables with structures
- Identification of optimal compounds or scaffolds
- Analysis of structure-activity trends
- Comparison of chemical series or scaffolds

**Typical user queries:**
- "Analyze the SAR of this compound library"
- "Identify scaffolds and perform R-group decomposition"
- "Which compounds have the best LLE?"
- "Generate SAR tables for each scaffold"
- "What are the key SAR trends in this dataset?"

---

## Required Input Data

The skill expects a **CSV file** with the following columns:

### Mandatory Columns:
- `Example_Number` or `Compound_ID`: Unique identifier for each compound
- `SMILES`: Molecular structure in SMILES notation
- Activity data (at least one):
  - `Ki_uM` or `IC50_uM` 

### Optional Columns:
- `Chemical_Name`: Compound name or identifier
- `cLogP`: Pre-calculated lipophilicity (will be calculated if missing)
- Any other experimental data

### Data Format Notes:
- Ki/IC50 values may contain special characters: `<`, `>`, `≥` (will be parsed correctly)
- Missing activity values are handled gracefully
- SMILES must be valid and parseable by RDKit

---

## Complete SAR Analysis Workflow

### Step 1: Initial Data Loading and Property Calculation

**Script**: `sar_analysis.py`

**Purpose**: Load the CSV data, parse SMILES structures, calculate molecular properties, and compute efficiency metrics.

**Usage**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py <input_csv>
```

**What it does**:
1. Loads CSV and validates SMILES structures
2. Parses Ki/IC50 values (handles <, >, ≥ symbols)
3. Calculates molecular properties:
   - Molecular weight (MW)
   - Lipophilicity (cLogP using Wildman-Crippen)
   - Hydrogen bond acceptors/donors (HBA/HBD)
   - Topological polar surface area (TPSA)
   - Rotatable bonds
   - Heavy atom count
4. Computes efficiency metrics:
   - pKi = -log₁₀(Ki in M)
   - Ligand Efficiency (LE) = pKi / Heavy Atoms
   - Lipophilic Ligand Efficiency (LLE) = pKi - cLogP
5. Generates summary statistics

**Output**: `data_with_analysis.csv` (full dataset with calculated properties)

**Example**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py Table1_Complete_with_SMILES.csv
```

---

### Step 2: Scaffold Identification and Clustering

**Script**: `scaffold_analysis.py`

**Purpose**: Generate Murcko scaffolds, cluster compounds by scaffold, and calculate scaffold-level statistics.

**Usage**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/scaffold_analysis.py
```

**Prerequisites**: Must run after `sar_analysis.py` (requires `data_with_analysis.csv`)

**What it does**:
1. Generates Murcko scaffolds (removes side chains, keeps ring systems)
2. Clusters compounds by identical scaffold
3. Identifies major scaffolds (≥5 compounds by default)
4. Calculates scaffold statistics:
   - Number of compounds per scaffold
   - Mean/median/min/max activity
   - Mean LLE and other properties
5. Assigns scaffold IDs (S1, S2, S3, etc.)

**Outputs**:
- `scaffold_statistics.csv`: Summary statistics per scaffold
- `scaffold_mapping.csv`: Compound-to-scaffold mapping
- `data_with_scaffolds.csv`: Full data with scaffold assignments

**Customization**: Edit the script to change the minimum compound threshold (default: 5)

---

### Step 3: R-Group Decomposition

**Script**: `rgroup_analysis.py`

**Purpose**: Decompose compounds into core scaffold + R-groups to understand substituent contributions.

**Usage**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/rgroup_analysis.py
```

**Prerequisites**: Must run after `scaffold_analysis.py`

**What it does**:
1. For each major scaffold, defines core structure with attachment points
2. Performs R-group decomposition using RDKit
3. Extracts R-group SMILES for each position
4. Maps R-groups to activity data
5. Identifies SAR trends per R-group position

**Output**: `rgroup_decomposition.csv`

**Columns**:
- `Core`: Core scaffold SMILES with attachment points ([*:1], [*:2], etc.)
- `R1`, `R2`, `R3`, ...: R-group SMILES for each position
- `Example_Number`: Compound identifier
- `Ki_BD1`, `Ki_BD2`, `Sum_Ki`: Activity data
- `cLogP`, `Avg_LLE`: Physicochemical properties
- `Scaffold_ID`: Scaffold identifier (S1, S2, etc.)

**Note**: The number of R-groups varies by scaffold (typically 4-9 positions)

---

### Step 4: Statistical Analysis and SAR Trends

**Script**: `sar_trends_analysis.py`

**Purpose**: Perform correlation analysis, identify top compounds, and analyze SAR trends.

**Usage**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_trends_analysis.py
```

**Prerequisites**: Requires `data_with_analysis.csv` and `scaffold_statistics.csv`

**What it does**:
1. Correlation analysis:
   - Pearson and Spearman correlations
   - cLogP vs potency, MW vs potency, etc.
2. Scaffold comparison:
   - One-way ANOVA or Kruskal-Wallis test
   - Identify significantly different scaffolds
3. Top compound identification:
   - Best potency (lowest Ki)
   - Best LLE (highest efficiency)
   - Best LE (highest ligand efficiency)
4. Generates summary statistics

**Output**: Console output with statistical results (can be redirected to file)

---

### Step 5: Visualization

**Script**: `create_visualizations.py`

**Purpose**: Generate publication-quality figures for SAR analysis.

**Usage**:
```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/create_visualizations.py
```

**Prerequisites**: Requires `data_with_analysis.csv` and `scaffold_statistics.csv`

**What it does**:
1. Potency vs lipophilicity scatter plot (color-coded by LLE)
2. Scaffold comparison bar charts (mean Ki, LLE)
3. Top compounds visualization
4. Property distribution histograms

**Outputs** (in `/figures/` directory):
- `potency_vs_lipophilicity.png`
- `scaffold_comparison.png`
- `top_potency_compounds.png`
- Additional figures as needed

**Customization**: Edit the script to adjust figure size, colors, or add additional plots

---

### Step 6: SAR Table Generation

**Scripts**: 
- `create_s1_s7_tables.py`: Generate detailed SAR tables for specific scaffolds
- `draw_all_scaffolds.py`: Generate scaffold structure images
- `draw_all_unique_scaffolds.py`: Comprehensive scaffold visualization

**Purpose**: Create publication-ready SAR tables with structures and activity data.

**Usage**:
```bash
# Generate SAR tables for major scaffolds
python3 /home/ubuntu/skills/sar-analyzer/scripts/create_s1_s7_tables.py

# Generate scaffold structure images
python3 /home/ubuntu/skills/sar-analyzer/scripts/draw_all_scaffolds.py

# Generate comprehensive scaffold overview
python3 /home/ubuntu/skills/sar-analyzer/scripts/draw_all_unique_scaffolds.py
```

**Prerequisites**: Requires `rgroup_decomposition.csv`

**What it does**:
1. Creates SAR tables with:
   - Core scaffold structure with R-group labels
   - R-group substituents for each compound
   - Activity data (Ki, LLE, cLogP)
   - Color-coded highlighting (best compounds in gold, excellent values in green)
   - Summary statistics
2. Generates high-resolution structure images
3. Creates scaffold comparison figures

**Outputs**:
- `SAR_Table_S1_Complete.png`, `SAR_Table_S2_Complete.png`, etc.
- `S1_Core_Structure.png`, `S2_Core_Structure.png`, etc.
- `all_scaffolds_batch*.png`
- `top20_scaffolds_summary.png`

---

## Efficiency Metrics Explained

### Ligand Efficiency (LE)

**Formula**: `LE = pKi / Heavy Atom Count`

**Interpretation**:
- LE > 0.3: Excellent
- LE 0.25-0.3: Good
- LE < 0.25: Poor

**Purpose**: Measures binding energy per atom; helps identify compact, efficient binders.

### Lipophilic Ligand Efficiency (LLE)

**Formula**: `LLE = pKi - cLogP`

**Calculation**:
1. Convert Ki (μM) to pKi: `pKi = -log₁₀(Ki × 10⁻⁶)`
2. Calculate LLE: `LLE = pKi - cLogP`

**Interpretation**:
- LLE > 7.0: Excellent (optimal drug-like properties)
- LLE 5.0-7.0: Good
- LLE 3.0-5.0: Moderate
- LLE < 3.0: Poor (high lipophilicity liability)

**Purpose**: Balances potency and lipophilicity; compounds with high LLE achieve good activity without excessive lipophilicity, reducing risk of poor PK and toxicity.

**Reference**: See `references/LLE_Calculation_Methodology.md` for detailed examples.

---

## Best Practices

### 1. Data Preparation
- **Validate SMILES**: Ensure all SMILES are valid before analysis
- **Standardize activity units**: Convert all Ki/IC50 to the same units (μM recommended)
- **Handle missing data**: Decide how to treat compounds with missing activity values
- **Check for duplicates**: Remove or flag duplicate structures

### 2. Scaffold Analysis
- **Adjust scaffold threshold**: For small datasets (<50 compounds), reduce minimum scaffold size to 3
- **Manual scaffold inspection**: Review scaffold structures to ensure they make chemical sense
- **Consider generic scaffolds**: Very generic scaffolds (e.g., single benzene ring) may not be informative

### 3. R-Group Decomposition
- **Define core carefully**: The core structure should represent the conserved scaffold
- **Attachment point numbering**: Ensure attachment points are numbered consistently
- **Handle complex R-groups**: Large or complex R-groups may need manual inspection
- **Missing R-groups**: Some compounds may not match the core perfectly; these will be excluded

### 4. Interpretation
- **LLE vs potency**: Prioritize compounds with high LLE (>7.0) over raw potency
- **Scaffold comparison**: Use statistical tests to confirm differences between scaffolds
- **SAR trends**: Look for consistent patterns across multiple R-group positions
- **Outliers**: Investigate compounds with unusual properties or activity

### 5. Visualization
- **Color coding**: Use consistent color schemes (e.g., green for good, red for poor)
- **Log scale**: Use log scale for activity data (Ki, IC50) in plots
- **Structure images**: Ensure structures are large enough to be readable
- **Publication quality**: Generate figures at 300 dpi for publication

---

## Common Workflows

### Workflow 1: Complete SAR Analysis (Recommended)

Run all scripts in sequence for comprehensive analysis:

```bash
# Step 1: Load data and calculate properties
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py input_data.csv

# Step 2: Identify scaffolds
python3 /home/ubuntu/skills/sar-analyzer/scripts/scaffold_analysis.py

# Step 3: R-group decomposition
python3 /home/ubuntu/skills/sar-analyzer/scripts/rgroup_analysis.py

# Step 4: Statistical analysis
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_trends_analysis.py

# Step 5: Generate visualizations
python3 /home/ubuntu/skills/sar-analyzer/scripts/create_visualizations.py

# Step 6: Create SAR tables
python3 /home/ubuntu/skills/sar-analyzer/scripts/create_s1_s7_tables.py
python3 /home/ubuntu/skills/sar-analyzer/scripts/draw_all_unique_scaffolds.py
```

### Workflow 2: Quick Property Calculation

For rapid property calculation without full SAR analysis:

```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py input_data.csv
```

Output: `data_with_analysis.csv` with all calculated properties

### Workflow 3: Scaffold-Focused Analysis

For projects focused on scaffold hopping or series comparison:

```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py input_data.csv
python3 /home/ubuntu/skills/sar-analyzer/scripts/scaffold_analysis.py
python3 /home/ubuntu/skills/sar-analyzer/scripts/draw_all_unique_scaffolds.py
```

### Workflow 4: R-Group SAR Tables Only

If scaffolds are already known and you want detailed R-group tables:

```bash
python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py input_data.csv
python3 /home/ubuntu/skills/sar-analyzer/scripts/scaffold_analysis.py
python3 /home/ubuntu/skills/sar-analyzer/scripts/rgroup_analysis.py
python3 /home/ubuntu/skills/sar-analyzer/scripts/create_s1_s7_tables.py
```

---

## Output Files Reference

### Data Files
| File | Description | Generated By |
|------|-------------|--------------|
| `data_with_analysis.csv` | Full dataset with calculated properties and metrics | `sar_analysis.py` |
| `scaffold_statistics.csv` | Scaffold-level summary statistics | `scaffold_analysis.py` |
| `scaffold_mapping.csv` | Compound-to-scaffold mapping | `scaffold_analysis.py` |
| `data_with_scaffolds.csv` | Full data with scaffold assignments | `scaffold_analysis.py` |
| `rgroup_decomposition.csv` | R-group analysis results | `rgroup_analysis.py` |

### Figure Files
| File | Description | Generated By |
|------|-------------|--------------|
| `potency_vs_lipophilicity.png` | Main SAR plot | `create_visualizations.py` |
| `scaffold_comparison.png` | Scaffold performance comparison | `create_visualizations.py` |
| `top_potency_compounds.png` | Best compounds visualization | `create_visualizations.py` |
| `SAR_Table_S*_Complete.png` | Detailed SAR tables per scaffold | `create_s1_s7_tables.py` |
| `S*_Core_Structure.png` | High-resolution scaffold images | `create_s1_s7_tables.py` |
| `all_scaffolds_batch*.png` | Comprehensive scaffold overview | `draw_all_unique_scaffolds.py` |
| `top20_scaffolds_summary.png` | Top 20 scaffolds summary | `draw_all_unique_scaffolds.py` |

---

## Troubleshooting

### Issue: RDKit cannot parse SMILES
**Solution**: Check for invalid SMILES strings in the input CSV. Use a SMILES validator or regenerate SMILES from structure files.

### Issue: No scaffolds identified
**Solution**: Reduce the minimum scaffold size threshold in `scaffold_analysis.py` (default: 5 compounds).

### Issue: R-group decomposition fails
**Solution**: Ensure the core structure is correctly defined with attachment points. Check that compounds actually match the core scaffold.

### Issue: Missing activity data
**Solution**: The scripts handle missing data gracefully. Compounds with missing Ki/IC50 will be excluded from specific analyses but retained in the dataset.

### Issue: Memory errors with large datasets
**Solution**: Process scaffolds separately or reduce image resolution in visualization scripts.

### Issue: Figures are too small/large
**Solution**: Edit the `figsize` parameter in visualization scripts (e.g., `figsize=(12, 8)`).

---

## Dependencies

The skill requires the following Python packages:

```
rdkit>=2023.9.1
pandas>=2.1.0
numpy>=1.24.3
matplotlib>=3.7.2
seaborn>=0.12.2
scipy>=1.11.2
```

**Installation**:
```bash
sudo pip3 install rdkit pandas numpy matplotlib seaborn scipy
```

**Note**: RDKit is pre-installed in the Manus sandbox environment.

---

## Advanced Usage

### Customizing Scaffold Definitions

To modify scaffold definitions or add custom cores:

1. Edit `rgroup_analysis.py`
2. Define custom core SMILES with attachment points:
   ```python
   core_smiles = "O=c1c2c(ccn1[*:5])C([*:1])=CC1=NN(c3ccc([*:3])cc3[*:4])C3C1=C2NCC3[*:2]"
   ```
3. Update the scaffold ID mapping

### Adding New Efficiency Metrics

To add custom efficiency metrics:

1. Edit `sar_analysis.py`
2. Add calculation after the LLE section:
   ```python
   df['Custom_Metric'] = df['pKi'] - (0.5 * df['cLogP']) + (0.1 * df['TPSA'])
   ```
3. Update visualization scripts to include the new metric

### Batch Processing Multiple Datasets

To analyze multiple datasets:

```bash
for file in *.csv; do
    echo "Processing $file"
    python3 /home/ubuntu/skills/sar-analyzer/scripts/sar_analysis.py "$file"
    # Add other scripts as needed
done
```

---

## References

For detailed methodology and additional information, see:

- `references/SAR_Analysis_Methods.md`: Complete methodology documentation
- `references/LLE_Calculation_Methodology.md`: Detailed LLE calculation guide
- `references/README_Scripts.md`: Quick reference for all scripts

---

## Example Use Case: BRD4 Inhibitor Analysis

**Scenario**: Analyze 276 BRD4 inhibitors to identify optimal scaffolds and R-groups.

**Input**: `Table1_Complete_with_SMILES.csv` with Ki_BDI_uM and Ki_BDII_uM data

**Steps**:
1. Run complete workflow (all 6 steps)
2. Review `scaffold_statistics.csv` to identify best scaffolds
3. Examine `rgroup_decomposition.csv` for R-group SAR trends
4. Generate SAR tables for top 3 scaffolds
5. Identify compounds with LLE > 7.0 for further testing

**Key Findings**:
- Scaffold S6 (pyrazole-fused): Best LLE (6.35), excellent potency
- Scaffold S3 (cyclohexyl): High LLE (6.37), consistent activity
- R1 position: CS(=O)(=O)C or CCS(=O)(=O)NH optimal
- Top compound: Example 22 (Sum Ki = 1.82 nM, LLE = 7.27)

---

## Skill Maintenance

This skill is designed to be modular and extensible. To add new functionality:

1. Add new Python scripts to `scripts/` directory
2. Update this `SKILL.md` with usage instructions
3. Add supporting documentation to `references/` if needed
4. Test with sample data before deployment

---

## Version History

- **v1.0.0** (2026-01-24): Initial release with complete SAR analysis workflow
