# Histone Code Figure Generation

This repository contains the scripts and local input files used to generate the current Figure 1, Figure 2, and Figure 3 panels.

The project is organized so that the main figure-generation scripts read from the local `Input/` folder structure and write outputs into `outputs/`.

## Input Layout

The current input layout is:

```text
Input/
├── Fig_1/
├── Fig_2/
└── Fig_3/
```

Additional Figure 1 metadata tables are kept at the project root:

- `Figure_1_A.csv`
- `Figure_1_B.csv`

## Main Generators

### Figure 1

Use:

```bash
python3 generate_figure1.py
```

Default inputs:

- `Figure_1_A.csv`
- `Figure_1_B.csv`
- `Input/Fig_1/`

Default output directory:

```text
outputs/figure1/
```

Generated outputs:

- `Figure_1_A_sequencing_availability.png`
- `Figure_1_A_legend.png`
- `Figure_1_B_states.png`
- `Figure_1_C.png`
- `Figure_1_D.png`
- `Figure_1_CD_legend.png`
- `Figure_1_E_E10_5.png`
- `Figure_1_E_P0.png`
- `Figure_1_E_legend.png`
- `Figure_1_F_emissions.png`
- `Figure_1_F_jaccard.png`

Useful CLI example:

```bash
python3 generate_figure1.py \
  --figure1-a-csv Figure_1_A.csv \
  --figure1-b-csv Figure_1_B.csv \
  --seg-input-dir Input/Fig_1 \
  --output-dir outputs/figure1 \
  --dpi 600 \
  --panel-width 3 \
  --panel-height 2
```

### Figure 2

Use:

```bash
python3 generate_figure2.py
```

Default inputs:

- `Input/Fig_1/` for Figure 2B emission summaries
- `Input/Fig_2/` for Figure 2C-F enrichment tables

Default output directory:

```text
outputs/figure2/
```

Useful CLI examples:

Generate all Figure 2 panels:

```bash
python3 generate_figure2.py
```

Generate only Figure 2B:

```bash
python3 generate_figure2.py --only-panel B
```

Generate one compact panel set:

```bash
python3 generate_figure2.py --compact-panel C
```

Custom input directories:

```bash
python3 generate_figure2.py \
  --fig1-input-dir Input/Fig_1 \
  --fig2-input-dir Input/Fig_2 \
  --output-dir outputs/figure2
```

### Figure 3

Use:

```bash
python3 generate_figure3.py
```

Default inputs:

- `Input/Fig_3/Adult_EO.txt`
- `Input/Fig_3/Ad_overlap.txt`
- `Input/Fig_3/EtoA_genome.txt`

Default output directory:

```text
outputs/figure3/
```

Generated outputs:

- `Fig3_AD_Emission.png`
- `Fig3_AD_FE.png`
- `Fig3D_gcover.png`

## Figure 3 Alternative Scripts

These scripts are kept as additional exploratory or replacement panels for Figure 3E:

- `Fig3E_transition_alternatives.py`
  Creates transition heatmap and bubble-matrix alternatives from:
  - `Input/Fig_3/transition_E10_E16.txt`
  - `Input/Fig_3/transition_E16_P0.txt`
  - `Input/Fig_3/transition_P0_adult.txt`

- `Fig3E_flow_alternatives.py`
  Creates simplified Sankey-style and alluvial-style flow plots from the same transition files.

- `Sankey_fig_3E.R`
  Original R Sankey workflow retained for reference.

Alternative outputs are written under:

- `outputs/figure3/transition_alternatives/`
- `outputs/figure3/transition_flows/`

## Standalone Panel Scripts

The repository also contains standalone plotting scripts for panel-level work and iterative editing, including:

- `plot_sequencing_timeline.py`
- `plot_figure_1_b.py`
- `plot_figure_1_cd.py`
- `plot_figure_1_ef.py`
- `plot_figure_1_f_original_style.py`
- `Fig_3_emission.py`
- `Fig3_FE.py`
- `Fig3D_gcover.py`

These are useful when adjusting one panel at a time, but the recommended entry points for reproducible generation are:

- `generate_figure1.py`
- `generate_figure2.py`
- `generate_figure3.py`

## Requirements

Python packages used across the scripts include:

```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn
```

Some legacy or optional workflows may also require:

```bash
pip install webshot2
```

The R Sankey script requires R packages such as:

- `dplyr`
- `tidyr`
- `readr`
- `jsonlite`
- `networkD3`
- `htmlwidgets`
- `webshot2`

## Recommended GitHub Upload Set

For a clean GitHub upload, the important reproducible pieces are:

- `generate_figure1.py`
- `generate_figure2.py`
- `generate_figure3.py`
- `Figure_1_A.csv`
- `Figure_1_B.csv`
- `Input/Fig_1/`
- `Input/Fig_2/`
- `Input/Fig_3/`

Optional extras:

- `Fig3E_transition_alternatives.py`
- `Fig3E_flow_alternatives.py`
- `Sankey_fig_3E.R`

## Notes

- Most scripts save at `600 dpi`.
- Figure 3 currently contains panel-specific styling refinements that were tuned directly in the current scripts.
- Some scripts may print non-fatal Matplotlib font-cache warnings depending on the local environment. These do not normally affect output generation.
