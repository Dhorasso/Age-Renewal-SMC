# Sequential Bayesian analysis of time-varying transmission in stratified epidemic renewal models

Code accompanying the paper *"Sequential Bayesian analysis of time-varying
transmission in stratified epidemic renewal models"*. The model is a
Bayesian age-stratified renewal model with time-varying, age-specific
susceptibility, fitted with a block SMC² algorithm.

## Repository structure

```
.
├── data/
│   ├── simulated/          # cached synthetic datasets (see load_or_simulate() below)
│   ├── real/                # Irish COVID-19 surveillance data (Covid_Data_Ireland.csv)
│   └── contact_matrices/    # age-structured contact matrices + age distribution (85-band CSVs)
├── src/                      # all model / inference code 
├── figures/
│   ├── real_data/            # figures generated from the real-data analysis
│   └── sim_data/             # figures generated from the simulation study
├── scripts/                   # entry-point scripts that call src/ on data/
└── README.md
```



## `scripts/` (entry points)

| Script | Fits | Data |
|---|---|---|
| `run_1age_covid.R` | 1-age EpiSSM | real total daily cases (`data/real/Covid_Data_Ireland.csv`) |
| `run_multi_simulation.R` | age-stratified EpiSSM | synthetic data, cached to `data/simulated/` |
| `run_ireland_analysis.R` | age-stratified EpiSSM | real age-stratified cases, with 14-day forecast + CRPS |
|`run_pmmh_ireland.R`      | age-stratified EpiSSM| inference via PMMH 

`run_pmmh_ireland.R` shows that PMMH() is an alternative to SMC2(): it reuses the exact same SSM and (almost) the exact
same opts as run_ireland_analysis.R, only adding the MCMC-specific fields (iterations, nChains, burnin, thin). The
same pattern works for any other model/data combination already set up for SMC2() in this repo — e.g. EpiSSM_1age on the real
1-age series, or EpiSSM on the simulated multi-age dataset: build opts the way the corresponding SMC2 script already does, add
the MCMC fields, and call PMMH(SSM, opts) instead of SMC2(SSM, opts). 
PMMH() requires the BayesianTools package (install.packages("BayesianTools")); it isn’t needed for any of the SMC2
based scripts.

Plot of the filtering-distribution are save in  `figures/real_data/` or `figures/sim_data/`.

**Load order.** in most scripts you only need:

```r
source("src/smc2.R")
source("src/posterior_marginal.R")
source("src/delay_distributions.R")
```

