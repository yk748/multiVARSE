# multi-VARSE

R code accompanying **“Joint modeling and inference of multiple-subject high-dimensional sparse vector autoregressive models”** by Younghoon Kim, Zachary F. Fisher, and Vladas Pipiras ([preprint, 2025](https://arxiv.org/abs/2510.14044)).

## Overview

This repository implements estimation and inference for high-dimensional vector autoregressive (VAR) models across multiple subjects. Each subject's transition matrix is decomposed into a **common component** shared across subjects and a **subject-specific (unique) component**.

The proposed framework:

- Fits sparse VAR models separately to individual subjects using Lasso.
- Debiases the individual estimates using nodewise regression.
- Identifies common and unique paths through robust aggregation based on a redescending loss, followed by thresholding.
- Selects aggregation and threshold parameters using cross-validation.
- Performs tests of **common-path significance**, **cross-subject homogeneity**, and **individual-path nullity**, with multiple-testing procedures.

The repository also includes comparison methods based on the stratified Lasso and adaptive stratified Lasso from the `multivar` R package, simulation experiments, and an fMRI data application.

## Repository structure

### Method implementation

| File | Purpose |
| --- | --- |
| `multiVARSE_optimizer.R` | Fits sparse subject-level VAR regressions using Lasso. |
| `multiVARSE_debiasor.R` | Computes debiased estimators via nodewise regression. |
| `multiVARSE_identifier.R` | Identifies common and unique components using aggregation and thresholding. |
| `multiVARSE_cv.R` | Cross-validates aggregation and threshold parameters. |
| `multiVARSE_tester.R` | Implements the three hypothesis tests and their multiple-testing procedures. |
| `multiVARSE_benchmark.R` | Runs the stratified Lasso and adaptive stratified Lasso benchmarks. |
| `multiVARSE_test_evaluator.R` | Evaluates false discovery rates and power in simulations. |

### Simulations and data analysis

| File | Purpose | Manuscript output |
| --- | --- | --- |
| `multiVARSE_simulator.R` | Generates sparse multi-subject VAR data. | Main simulations and higher-dimensional simulations |
| `multiVARSE_simulator_overlap.R` | Generates data under alternative overlap configurations. | Supplementary simulations |
| `multiVARSE_parallel_main.R` | Runs the main parallel simulations, summarizes estimation/testing performance, and produces plots. | Figures 1, 2, and S.1 |
| `multiVARSE_parallel_higher_dim.R` | Runs higher-dimensional simulation experiments. | Figure S.2 |
| `multiVARSE_parallel_additional.R` | Runs the additional overlap simulations. | Figure S.3 |
| `multiVARSE_data_application.R` | Fits the benchmark and proposed methods to processed fMRI time series and visualizes estimated paths. | Figure 3 |

## Requirements

The code is written in R. Install the packages used by the simulation, estimation, visualization, and data-application scripts:

```r
install.packages(c(
  "mvtnorm", "glmnet", "Matrix", "multivar", "magic",
  "doParallel", "foreach", "dqrng", "doRNG",
  "ggplot2", "latex2exp", "patchwork", "RColorBrewer", "reshape2"
))
```

The optional preprocessing code for raw fMRI files additionally uses `ciftiTools` and Connectome Workbench. It is commented out in the supplied data-application script.

## Reproducing the results

Run the scripts from the repository's root directory so their relative `source()`, `load()`, and data-file paths resolve correctly.

### Main simulation results

```bash
Rscript multiVARSE_parallel_main.R
```

This script uses a VAR(1) model across the following settings:

- Number of subjects: `K = 10, 15`
- Dimension: `d = 10, 20`
- Mean time-series length: `100, 200`
- Sparsity scenarios: `case1`, `case2`, `case3`
- Monte Carlo iterations: `50` per configuration

The script saves simulation results as files with names such as `multivarse_d10_T100_subject10_case1.RData` and generates `figure1.pdf`, `figure2.pdf`, and `figureS1.pdf`.

### Additional simulations

```bash
Rscript multiVARSE_parallel_higher_dim.R
Rscript multiVARSE_parallel_additional.R
```

- **Higher-dimensional simulations:** `d = 30, 50`, `K = 10`, mean series length `100` or `200`, and three sparsity scenarios; produces `figureS2.pdf`.
- **Additional overlap simulations:** `d = 10, 20`, `K = 10`, mean series length `100` or `200`, and three overlap scenarios; produces `figureS3.pdf`.

Both scripts use 50 Monte Carlo iterations per configuration by default and save intermediate `.RData` files. The simulation scripts create parallel workers, using up to 50 cores while reserving one detected core when possible; adjust the worker settings and number of iterations for your machine before running large experiments.

### fMRI data application

```bash
Rscript multiVARSE_data_application.R
```

The script expects **12 preprocessed subject-level CSV files** named `HCP_Emotion_<subject_id>.csv` in the working directory. Each file contains time-series measurements for 17 functional networks. The script loads these files, estimates the path structures, performs the three tests, and generates `figure3.pdf`.

The script includes a commented preprocessing workflow for deriving network-level time series from CIFTI files. To reproduce the application, first ensure that the required processed CSV files are available at the expected paths. Raw imaging data are not needed when using these processed files.

## Notes on execution

- The parallel simulation scripts execute their full simulation grids when sourced or run, so they may take substantial time and memory.
- They write `.RData` and PDF outputs into the current working directory and can overwrite files with the same names.
- Parallel iterations that fail are omitted by the current `foreach` configuration (`.errorhandling = "remove"`); check the reported number of successful iterations before interpreting results.
- The code has been reviewed for its documented structure, but the full simulations and fMRI analysis have **not** been executed as part of this README update.

## References

1. Kim, Y., Fisher, Z. F., & Pipiras, V. (2025). *Joint modeling and inference of multiple-subject high-dimensional sparse vector autoregressive models*. [arXiv:2510.14044](https://arxiv.org/abs/2510.14044).
2. Fisher, Z. F., Kim, Y., Pipiras, V., Crawford, C., Petrie, D. J., Hunter, M. D., & Geier, C. F. (2024). Structured estimation of heterogeneous time series. *Multivariate Behavioral Research, 59*(6), 1270–1289. [DOI](https://doi.org/10.1080/00273171.2023.2283837).
3. Fisher, Z. F., Kim, Y., Fredrickson, B. L., & Pipiras, V. (2022). Penalized estimation and forecasting of multiple subject intensive longitudinal data. *Psychometrika, 87*(2), 403–431. [Article](https://www.cambridge.org/core/journals/psychometrika/article/abs/penalized-estimation-and-forecasting-of-multiple-subject-intensive-longitudinal-data/A5933DB71583B091EB3D65CBA2E28498).
