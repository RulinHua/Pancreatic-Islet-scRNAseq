# Pancreatic-Islet-scRNAseq

**Author:** Rulin Hua  
**Primary language:** R, with Bash/Python for pySCENIC

Reproducible computational workflows for human pancreatic islet single-cell
RNA-seq analysis, with an emphasis on disease/HbA1c stratification and
donor-level beta-cell analyses.

> **Manuscript status:** this repository contains analysis code associated
> with ongoing, unpublished work. Analysis outputs and biological
> interpretation may change as the study develops.

## Analysis modules

- Cell-type UMAPs, composition summaries, and marker visualization
- SCN-family expression across Control/T1D/T2D and Normal/Prediabetes/Diabetes
- Donor-level pseudobulk differential-expression analysis with DESeq2
- Covariate-aware SCN-family/TRPC4 pseudobulk models with focused beta-cell SCN9A and alpha-cell TRPC4 analyses
- Continuous donor-level SCN9A association models in beta cells
- Top-gene/module heatmaps and GO/KEGG GSEA
- Predefined calcium/exocytosis functional-module analysis
- CytoTRACE2 beta-cell potency analysis in HPAP with IIDP replication
- pySCENIC regulatory-network inference and donor-level regulon analysis
- SCN9A-associated upstream TF/regulon network analysis
- HbA1c-informed Monocle2 beta-cell pseudotime analysis

## Repository layout

```text
.
├── R/
│   └── helpers.R
├── config/
│   └── paths.example.R
├── data/
│   └── README.md
├── results/
│   └── README.md
├── scripts/
│   ├── 00_setup.R
│   ├── 01_celltype_overview.R
│   ├── 02_SCN_family_expression.R
│   ├── 03_pseudobulk_DESeq2.R
│   ├── 04_targeted_SCN_pseudobulk.R
│   ├── 05_SCN9A_beta_cell_association.R
│   ├── 06_GSEA_and_functional_modules.R
│   ├── 07_CytoTRACE2.R
│   ├── 08_SCENIC_prepare.R
│   ├── 09_run_pyscenic.sh
│   ├── 10_SCENIC_downstream.R
│   ├── python/
│   │   └── 11_export_regulon_targets.py
│   ├── 12_SCENIC_networks_and_exploratory_monocle3.R
│   └── 13_monocle2_pseudotime.R
├── .gitignore
└── RUN_ORDER.md
```

## Input data

The single-cell RNA-seq data used in this analysis were obtained from
[PanKbase](https://data.pankbase.org/).

The analysis is based on the PanKbase resource analysis set
[PKBDS1349YHGQ](https://data.pankbase.org/analysis-sets/PKBDS1349YHGQ/),
which contains a reference map of human pancreatic islet cell-type-specific
gene expression derived from HPAP, IIDP, and Prodo samples.

The corresponding PanKbase matrix file is
[PKBFI5903OGWY](https://data.pankbase.org/matrix-files/PKBFI5903OGWY/).

For the analyses in this repository, the data were loaded from a local Seurat
RDS object (`060425_scRNA_v3.3.rds`). The large data object is not redistributed
in this repository.

Copy `config/paths.example.R` to `config/local_paths.R` and edit the local
input/output paths. The local configuration file is excluded by `.gitignore`.

## Running the workflow

This analysis was developed interactively rather than as a single monolithic
pipeline. Run:

```r
source("scripts/00_setup.R")
source("scripts/01_celltype_overview.R")
```

and continue through the numbered scripts in the same R session. See
`RUN_ORDER.md` for the complete sequence.

pySCENIC is run separately from a conda environment:

```bash
conda activate SCENIC
export MOTIF_DIR=/path/to/cisTarget_databases_hg38
bash scripts/09_run_pyscenic.sh
python scripts/python/11_export_regulon_targets.py
```

## Software

Major R/Bioconductor packages used include Seurat, data.table, DESeq2,
clusterProfiler, fgsea, GSVA, limma, ComplexHeatmap, CytoTRACE2, Monocle,
igraph/ggraph, and supporting visualization packages.

The SCENIC workflow additionally requires pySCENIC, loom support, the human
transcription-factor list, motif annotations, and hg38 cisTarget databases.

## Reproducibility notes

- Donor-level pseudobulk counts are aggregated from raw single-cell counts.
- Disease/HbA1c comparisons require minimum donor representation.
- The SCENIC GRN step balances HPAP donors by subsampling cells per donor.
- IIDP V3 beta cells are used as an independent validation cohort in several
  downstream analyses.
- Dataset-specific metadata harmonization from the original workflow is
  retained in `scripts/00_setup.R`.

## Data and outputs

Raw/processed single-cell objects and generated result files are excluded
from git. This keeps the repository focused on analysis code and avoids
redistributing large third-party datasets.

## License / reuse

No open-source license is granted at this stage because the associated
manuscript is ongoing. The code is publicly viewable for transparency and
portfolio/reproducibility purposes. Please contact the author before reuse.
