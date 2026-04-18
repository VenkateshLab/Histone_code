# Histone_code Figure Package

This folder is arranged for GitHub upload with only two runnable scripts and the required inputs.

## Folder structure

- `code/generate_figure1.py`: standalone script for all Figure 1 panels
- `code/generate_figure2.py`: standalone script for all Figure 2 panels
- `input/Figure_1_A.csv`: Figure 1A input table
- `input/Figure_1_B.csv`: Figure 1B input table
- `input/Fig_1/`: Figure 1 ChromHMM segment and emission files
- `input/Fig_2/`: Figure 2 fold-enrichment and enrichment files
- `outputs/`: default output location when the scripts are run

## Usage

Run Figure 1:

```bash
python3 code/generate_figure1.py
```

Run Figure 2:

```bash
python3 code/generate_figure2.py
```

## Notes

- Figure 1 outputs are written to `outputs/figure1/`
- Figure 2 outputs are written to `outputs/figure2/`
- Both scripts use the local `input/` folder by default, so the package can be moved as one folder
