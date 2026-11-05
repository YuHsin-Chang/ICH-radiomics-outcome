# ICH radiomics outcome prediction — data and reproduction code

Supporting data and analysis code for:

> Chang Y-H, Chen P-C, et al. An automated CT-based radiomics pipeline for functional outcome prediction in intracerebral hematoma: a feasibility study. *BMC Medical Imaging* (under review).

## Contents

| File | Description |
|---|---|
| `ICH_radiomics_deidentified.xlsx` | De-identified data for the 67 patients in the outcome analysis. Sheets: `clinical` (outcome and clinical predictors), `ground_truth` (107 PyRadiomics features from the adjudicated physician-defined segmentation), `doctor1`–`doctor3` (features from each neurosurgeon's independent segmentation, used for ICC), `DL_predict` (features from the deep-learning segmentation). Patient IDs (P001–P067) are study codes with no link to clinical identifiers. |
| `reproduce_models.R` | Reproduces feature-level ICC filtering, correlation filtering, LASSO feature selection, Rad-score construction, the hierarchical logistic models (Models 1–3), DeLong comparisons, and bootstrap optimism correction. |

### `clinical` sheet

| Column | Definition |
|---|---|
| `ID` | Study code |
| `outcome_poor_mRS4to6` | 1 = modified Rankin Scale 4–6 at hospital discharge; 0 = mRS 0–3 |
| `age` | Age in years |
| `sex_male` | 1 = male; 0 = female |
| `BMI` | Body mass index, kg/m² |
| `GCS` | Glasgow Coma Scale score at initial emergency department evaluation |
| `DL_volume_mL` | Hematoma volume (mL) from the deep-learning segmentation |

Radiomic features were extracted with PyRadiomics 3.0.1 from the original Hounsfield-unit images (1 × 1 × 1 mm resampling, fixed bin width 5 HU, no normalization, no filters); see the Supplementary Methods of the manuscript.

## Running the code

```r
install.packages(c("readxl", "irr", "glmnet", "pROC"))
setwd("<folder containing both files>")
source("reproduce_models.R")
```

Tested with R 4.5.1, glmnet 4.1-10, pROC 1.19.0.1. The script prints Table 3, the DeLong p-values, the Model 3 coefficients (Supplementary Table S4), and the bootstrap results (Supplementary Tables S5–S6). Bootstrap estimates vary slightly between runs because of resampling; with the default seed the script gives an optimism-corrected AUC of 0.838 for Model 3 (manuscript: 0.841) and feature selection frequencies within ±5 percentage points of the reported values.

## Not included

CT images, segmentation masks and the trained segmentation model cannot be shared because of institutional restrictions. An executable version of the segmentation model can be made available to collaborators under a collaboration agreement; contact the corresponding author.

## License

Data: CC BY 4.0. Code: MIT.
