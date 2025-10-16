# nature_in_vivo_ht7

This repository contains two analysis components used in our Nature manuscript:
1) **gwas_ko/**: pooled CRISPR knockout screen across multiple donors, analyzed with MAGeCK.
2) **tcr_analysis/**: TCR repertoire analysis from MiXCR outputs.

## Quickstart

### Environment
```bash
mamba env create -f environment.yml
mamba activate nature_in_vivo_ht7
```

### GWAS / KO screen
- Place FASTQs and update `gwas_ko/data/sample_info.csv` with columns: `sample_id, donor, bin, r1, r2`.
- Run (dry run to print commands):
```bash
cd gwas_ko
python scripts/run_mageck_clean.py --sample-info data/sample_info.csv --library data/brunello_library.txt --outdir ../results/mageck --dry-run
```
- Open notebook:
```
jupyter lab gwas_ko/notebooks/MAGeCK_GW_analysis-4_clean.ipynb
```

### TCR analysis
- Follow docs in `tcr_analysis/docs/MIXCR_QUICKSTART.md` to produce MiXCR outputs.
- Open notebook:
```
jupyter lab tcr_analysis/notebooks/TCR_090925_clean.ipynb
```

## Data policy
Large raw data are not stored in this repository. Only small reference files are included for reproducibility. See manuscript methods for data access.

## Citation
See `CITATION.cff` and the Zenodo DOI once archived.