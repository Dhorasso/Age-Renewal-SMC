# Sequential Bayesian analysis of time-varying transmission in stratified epidemic renewal models

Code accompanying the paper *"Sequential Bayesian analysis of time-varying transmission in stratified epidemic renewal models"*. The model is a Bayesian age-stratified renewal model with time-varying, age-specific susceptibility, fitted with a block SMC² algorithm.

## Getting started

### 1. Clone the repository

Clone the GitHub repository to your local machine using Git:

```bash
git clone https://github.com/Dhorasso/Age-Renewal-SMC.git
```

Navigate into the project directory:

```bash
cd Age-Renewal-SMC
```

Alternatively, you can download the repository as a ZIP file from the [GitHub repository](https://github.com/Dhorasso/Age-Renewal-SMC) by clicking **Code → Download ZIP**.

### 2. Repository structure

```text
.
├── data/
│   ├── simulated/          # Cached synthetic datasets
│   ├── real/               # Irish COVID-19 surveillance data
│   └── contact_matrices/   # Contact matrices and age distribution
├── figures/
│   ├── real_data/          # Figures from the real-data analysis
│   └── sim_data/           # Figures from the simulation study
├── scripts/                # Entry-point scripts that call src/ on data/
├── src/                    # Model and inference code
├── LICENSE
└── README.md
```

## `scripts/` (entry points)

| Script | Model fitted | Data |
|---|---|---|
| `run_multi_simulation.R` | Age-stratified EpiSSM | Synthetic data, cached to `data/simulated/` |
| `run_ireland_analysis.R` | Age-stratified EpiSSM | Real age-stratified cases, with a 14-day forecast |
| `run_pmmh_ireland.R` | Age-stratified EpiSSM | Inference via PMMH |
| `run_1age_covid.R` | Homogeneous EpiSSM | Real total daily cases (`data/real/Covid_Data_Ireland.csv`) |
| `run_epinow2.R` | EpiNow2 model | Real total daily cases (`data/real/Covid_Data_Ireland.csv`) |

## Alternative inference using PMMH

`run_pmmh_ireland.R` demonstrates that `PMMH()` can be used as an alternative to `SMC2()`.

It reuses the same state-space model (SSM) and almost the same options as `run_ireland_analysis.R`, adding only the MCMC-specific settings: `iterations`, `nChains`, `burnin`, and `thin`.

The same approach can be applied to other model/data combinations already configured for `SMC2()` in this repository. For example:

- `EpiSSM_1age` on the real one-age-group series.
- `EpiSSM` on the simulated multi-age dataset.

To switch from SMC² to PMMH, construct the options as in the corresponding SMC² script, add the MCMC-specific fields, and call `PMMH(SSM, opts)` instead of `SMC2(SSM, opts)`.

PMMH requires the `BayesianTools` R package:

```r
install.packages("BayesianTools")
```

This package is not required for scripts that use SMC² only.

## Output figures

Plots of the filtering distributions are saved in:

- `figures/real_data/`
- `figures/sim_data/`

## Loading the source files

Most scripts require the following source files:

```r
source("src/smc2.R")                # Or source("src/pmmh.R")
source("src/marginal_state.R")      # Or source("src/posterior_marginal.R")
source("src/delay_distributions.R")
```

Run the commands from the root directory of the cloned repository so that the relative file paths resolve correctly.
