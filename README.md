# ROSMAP Multi-Omic DIABLO Analysis
### Integrating proteomics and metabolomics to identify Alzheimer's disease biomarkers

---

## What this is

This repository contains the analysis pipeline for a narrative review examining multi-omic signatures across the cognitive continuum in Alzheimer's disease — from no cognitive impairment (NCI) through mild cognitive impairment (MCI) to AD. The analysis uses the ROSMAP cohort and integrates brain proteomics with brain metabolomics using the DIABLO framework from the mixOmics R package.

The goal was not a primary study but a structured integration of existing ROSMAP data to surface cross-omic biomarker patterns that are interpretable in the context of the existing AD literature.

---

## Dataset

**ROSMAP (Religious Orders Study and Memory and Aging Project)**
Data accessed via the AD Knowledge Portal (Synapse). You need a Synapse account and approved data access to download the raw files — they are not included here for data privacy reasons.

| Block | File used | Samples |
|---|---|---|
| Proteomics | `ROSMAP_proteomics_merged_clean.csv` | 340 |
| Metabolomics | `ROSMAP_Metabolon_HD4_Brain514_assay_data.csv` | 521 |
| Clinical | `ROSMAP_clinical.csv` | 3584 |
| **Matched overlap** | — | **225** |

**Outcome variable:** `cogdx` — consensus cognitive diagnosis
- 1 = NCI (n=95), 2 = MCI (n=59), 4 = AD (n=71)
- cogdx groups 3, 5, 6 excluded (insufficient sample size)

> Gene expression and genetic variant blocks were not included. The sample ID mapping files needed to link those datasets were not available in this data export.

---

## Method

DIABLO (`block.splsda`) from the mixOmics package. DIABLO is a supervised multi-block integration method that simultaneously selects features across omic layers that discriminate between groups. It returns sparse canonical variates — essentially, the smallest set of proteins and metabolites that together best separate NCI, MCI, and AD.

Key settings:
- Design matrix: 0.1 (weak inter-block correlation prior)
- Components: 2
- Feature selection tuned via 5-fold cross-validation, 10 repeats
- Circos plot correlation cutoff: r = 0.6

---

## How to reproduce this

**1. Get the data**
Download the files listed above from Synapse and place them in the folder structure shown below.

**2. Install R dependencies**
```r
install.packages(c("readr", "dplyr", "BiocManager"))
BiocManager::install("mixOmics")
```

**3. Update the file paths**
Each script has a `base_var`, `base_prot`, and `base_met` variable at the top. Point these at wherever you downloaded the data.

**4. Run the scripts in order**
```r
source("scripts/01_load_and_match.R")
source("scripts/02_preprocess.R")
source("scripts/03_diablo_model.R")   # takes 5-15 min
source("scripts/04_figures.R")
source("scripts/05_session_info.R")
```

---

## Folder structure

```
rosmap_diablo/
├── README.md
├── .gitignore
├── scripts/
│   ├── 01_load_and_match.R       # sample ID alignment across blocks
│   ├── 02_preprocess.R           # matrix building, imputation, filtering
│   ├── 03_diablo_model.R         # tuning, model fitting, performance
│   ├── 04_figures.R              # Figure 1 (sample plots) + Figure 2 (circos)
│   └── 05_session_info.R         # package versions used
├── results/
│   ├── 03_diablo_model.RData     # fitted model object
│   ├── biomarker_summary.csv     # selected features with loadings
│   └── session_info.txt          # R + package versions
└── figures/
    ├── Figure1_sample_plots.png
    └── Figure2_circos.png
```

---

## Key findings

The DIABLO model identified a cross-omic signature centred on three converging themes:

1. **Membrane phospholipid dysregulation** — elevated glycerophospholipids in AD reflecting accelerated membrane breakdown
2. **Antioxidant depletion** — progressive loss of gamma-tocopherol and homocarnosine across NCI → MCI → AD
3. **Synaptic and axonal protein loss** — NRN1, SNAP25, and PLXNB1 all decline toward AD, correlated with the metabolic shifts above

MCI shows an intermediate metabolic profile for most features, though three markers (N-acetyl-3-methylhistidine, homocarnosine, gamma-tocopherol) decline gradually and consistently across all three groups, suggesting potential utility as early trajectory markers.

---

## Limitations

- 225 samples after cross-block matching — modest for a 3-class problem
- MCI classification performance was poor (error rate ~90%), consistent with known clinical heterogeneity of MCI
- Analysis is cross-sectional; longitudinal validation would strengthen the biomarker candidates
- Gene expression and variant data could not be integrated due to missing ID mapping files

---

## R package versions

See `results/session_info.txt` for the full environment snapshot.
Core packages: mixOmics 6.34.0, R 4.x

---

## Citation

If you use this pipeline, please cite the mixOmics package:

> Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: An R package for 'omics feature selection and multiple data integration. *PLOS Computational Biology*, 2017.

And the ROSMAP study:

> Bennett DA et al. Religious Orders Study and Rush Memory and Aging Project. *Journal of Alzheimer's Disease*, 2018.
