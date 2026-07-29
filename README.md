# Heat-Health Warning System and Hospital Admissions in Germany

This repository contains the R code and technical documentation used for the analyses conducted as part of the master's thesis evaluating the German Heat-Health Warning System and hospital admissions.

The workflow covers:

- preparation of district-level environmental data;
- reconstruction of pre-implementation heat-alert classifications;
- preparation of the input dataset for analyses in the secure Research Data Center (FDZ) environment;
- district-specific first-stage analyses;
- second-stage meta-analyses;
- sensitivity and robustness analyses; and
- generation of tables and figures.

All scripts use project-relative paths and are intended to be run from the repository root.

## Repository contents

```text
project-root/
├── README.md
├── DWD_Warning_Cell_to_Harmonized_AGS_Crosswalk.xlsx
├── ICD10_Outcome_Variable_Dictionary.pdf
├── r/
│   ├── 01_exposure_preparation/
│   ├── 02_data_assembly_and_predictor_selection/
│   ├── 03_random_forest_models/
│   ├── 04_fdz_input_preparation/
│   ├── 05_fdz_analyses/
│   └── 06_fdz_output_processing/
├── data/
│   ├── raw/
│   ├── intermediate/
│   ├── processed/
│   ├── external/
│   └── restricted/
└── outputs/
    ├── random_forest_models/
    ├── 05_fdz_analyses/
    └── 06_fdz_output_processing/
```

Restricted data, large raw input files, and disclosure-controlled FDZ outputs are not included in the public repository.

## Software

All analyses were conducted in R version 4.3.1 (R Core Team, 2023) using RStudio version 2023.09.1+494 (Posit Team, 2023).

Required R packages are listed in the header of each script.

### Software references

- R Core Team. (2023). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>
- Posit Team. (2023). *RStudio: Integrated Development Environment for R*. Posit Software, PBC. <https://posit.co/>

## Required data

The workflow requires:

- gridded meteorological data;
- harmonized German district boundaries;
- 2011 Census population data;
- official DWD heat-warning data;
- a state-level public-holiday table for 2000–2022;
- restricted annual hospital-admission data; and
- disclosure-checked FDZ result files.

The main data locations are:

```text
data/raw/                         Raw environmental, spatial, population, and calendar data
data/intermediate/                Derived exposure and covariate files
data/processed/random_forest_input/
data/external/fdz_input/          Prepared FDZ input files
data/restricted/                  Restricted hospital data
data/processed/fdz_outputs/       Released FDZ result tables
```

The exact filenames, required variables, inputs, and outputs are documented in the header of each script.

## Recommended execution order

### 1. Environmental exposure preparation

Run all scripts in:

```text
r/01_exposure_preparation/
```

in numerical order.

These scripts prepare daily district-level meteorological variables, population-weighted exposures, and apparent temperature.

### 2. Data assembly and predictor selection

Run:

```text
02_01_assemble_annual_climate_exposure_data.R
02_02_prepare_random_forest_input_dataset.R
02_03_assess_predictor_correlations.R
```

`02_03` uses the object created by `02_02` and should therefore be run in the same R session unless it is adapted to load the saved dataset directly.

### 3. Random Forest models

Run:

```text
03_01_fit_over_and_undersampled_random_forest_models.R
03_02_fit_random_forest_without_population_weighted_predictors.R
03_03_fit_weighted_random_forest_and_evaluate_thresholds.R
```

The resulting predictions for 2000–2004 are used to reconstruct pre-implementation alert-qualifying days.

The primary alert definition used in the main analyses is `heat_alerts_rfm5`, based on the weighted Random Forest model with a probability threshold of 0.80. Other model variants are retained for sensitivity analyses.

### 4. FDZ input preparation

Run:

```text
04_create_district_holiday_calendar.R
04_create_fdz_exposure_and_alert_master_dataset.R
```

The first script creates a daily holiday indicator for each district. The second combines:

- meteorological exposures;
- holiday indicators;
- official DWD warnings; and
- reconstructed pre-implementation alert classifications.

Expected outputs:

```text
data/external/fdz_input/Masterfile_Final_for_FDZ.rds
data/external/fdz_input/Masterfile_Final_for_FDZ.csv
```

### 5. Secure FDZ analyses

The scripts in `r/05_fdz_analyses/` are organized as configuration-analysis pairs:

```text
05_01 → 05_02   Create district-specific daily time series
05_03 → 05_04   All-cause first-stage models
05_05 → 05_06   Cause-specific first-stage models
05_07 → 05_08   Separate pre- and post-period models
05_09 → 05_10   Post-period reference models
05_11 → 05_12   Combined first-stage and DF-sensitivity models
```

Each configuration script defines paths, model parameters, outcomes, analysis periods, logging settings, and the paired analysis script.

Before execution, verify:

```r
syntax_path <- "r/05_fdz_analyses"
```

and ensure that `syntax_name` exactly matches the paired script filename without the `.R` extension.

Restricted hospital data are expected under either:

```text
data/restricted/local/
```

or:

```text
data/restricted/fdz/
```

depending on the selected environment.

### 6. FDZ output processing

Place disclosure-checked FDZ output files in:

```text
data/processed/fdz_outputs/
```

Then run the main scripts in:

```text
r/06_fdz_output_processing/
```

Sensitivity and robustness analyses are located in:

```text
r/06_fdz_output_processing/06_07_sensitivity_and_robustness/
```

Additional descriptive and exploratory figures are created with:

```text
06_08_create_additional_descriptive_and_robustness_figures.R
```

`06_08` uses the district list produced by `06_02`, so `06_02` must be run first when that subgroup analysis is required.

## Reproducibility notes

- Run all scripts from the repository root.
- Some scripts create their output directories automatically, whereas others assume that the required directories already exist.
- Files in `data/processed/fdz_outputs/` must retain the filenames referenced in the scripts unless all corresponding paths are updated.
- Some scripts represent alternative or sensitivity analyses and are not required for reproduction of the primary estimates.
- Filenames, paths, script headers, and comments were standardized for repository publication without intentionally modifying the statistical analyses.

## Technical reference files

Three technical reference files accompany the code:

### `DWD_Warning_Cell_to_Harmonized_AGS_Crosswalk.xlsx`

This workbook contains the complete mapping of German Meteorological Service warning-cell identifiers to the harmonized 2022 administrative district codes used in the analysis.

It documents:

- direct correspondences;
- administrative reforms;
- aggregation rules;
- relevant years; and
- rationales for manually reviewed assignments.

The crosswalk was used to assign official DWD heat warnings to the harmonized districts used in the epidemiological analyses.

### `ICD10_Outcome_Variable_Dictionary.pdf`

This file documents the ICD-10-derived hospital-admission outcome variables generated during data preparation.

It includes:

- variable names;
- outcome definitions;
- included ICD-10 codes;
- excluded ICD-10 codes; and
- whether each variable was included in the final analyses.

### `public_holidays_2000-2022.csv`

This file contains all nationwide and state-specific statutory public holidays in Germany from 2000 through 2022, including the date, holiday name, and applicable federal state or region. The file was used to construct the district-day public-holiday indicator included in the time-series models.

## Data availability

The restricted hospital data are not publicly available and can be accessed only through the Research Data Centers of the Federal Statistical Office and the statistical offices of the Fedeal States under the applicable access and disclosure-control procedures.

Meteorological data, official heat-warning data, administrative boundary data, and population data were obtained from the publicly accessible official sources cited in the associated thesis or publication. Access and reuse remain subject to the respective providers' terms of use.

The analysis code and documentation of the data-processing workflow are provided in this repository. Restricted hospital data, large raw input files, and disclosure-controlled FDZ outputs are not included.
