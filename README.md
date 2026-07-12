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
├── src/                      # all model / inference code (see module map below)
├── figures/
│   ├── real_data/            # figures generated from the real-data analysis
│   └── sim_data/             # figures generated from the simulation study
├── scripts/                   # entry-point scripts that call src/ on data/
└── README.md
```


Only the overall contact matrix (`Ireland_country_level_M_overall_contact_matrix_85.csv`)
and the age distribution are currently used by the scripts. The four
setting-specific matrices (community/household/school/work) are stored
in `data/contact_matrices/` for future use but aren't wired into any
driver yet.

## `src/` module map

| File                      | Contents |
|---------------------------|----------|
| `utils.R`                 |  `logSumExp`, progress bar, `par_lapply` |
| `resampling.R`             | `Resample()`: systematic / multinomial / stratified |
| `particle_filter.R`       | Bootstrap particle filter: `OneStepPrediction`, `OneStepUpdate`, `BootstrapPF` |
| `delay_distributions.R`   | Discretized reporting-delay / generation-time PMFs |
| `priors.R`                | Prior sampling and (log-)density evaluation |
| `mcmc_kernels.R`          | `HWG_kernel`, `DA_HWG_kernel`, kNN surrogate, adaptive `Nx` |
| `smc2.R`                  | `IncrementalPF`, `PF_SMC2`, main `SMC2()` loop |
| `posterior_marginal.R`   | `PosteriorMarginal()`: posterior sample of latent trajectories |
| `epi_ssm_1age.R`          | Single-population EpiSSM (`InitState`, `StateProcess`, `ObsProcess`) |
| `epi_ssm_multiage.R`      | Age-stratified EpiSSM, contact-matrix-coupled renewal model |
| `contact_matrix_utils.R`  | Reciprocity correction + age-group aggregation of contact matrices |
| `epi_diagnostics.R`       | Next-generation-matrix R_t / R_{a,t}, aggregate population |
| `diagnostics.R`           | All ggplot2 plotting helpers (filtering ribbons, SMC² diagnostics, prior/posterior overlays, grid layout) |
| `sim_data_1age.R`         | Generation-time / report-delay kernels for the 1-age driver |
| `sim_data_multi.R`        | Forward-simulates synthetic age-stratified data for the simulation study |
| `simulated_data_cache.R`  | `load_or_simulate()`: cache any simulator's output to CSV under `data/simulated/` |
|`pmmh.R`                    | PMMH(): Particle Marginal Metropolis-Hastings via BayesianTools package, using the same `BootstrapPF` likelihood|

## `scripts/` (entry points)

| Script | Fits | Data |
|---|---|---|
| `run_1age_covid.R` | 1-age EpiSSM | real total daily cases (`data/real/Covid_Data_Ireland.csv`) |
| `run_multi_simulation.R` | age-stratified EpiSSM | synthetic data, cached to `data/simulated/` |
| `run_ireland_analysis.R` | age-stratified EpiSSM | real age-stratified cases, with 14-day forecast + CRPS |
|`run_pmmh_ireland.R`      | age-stratified EpiSSM, via PMMH |

run_pmmh_ireland.R shows that PMMH() is an alternative to SMC2(): it reuses the exact same SSM and (almost) the exact
same opts as run_ireland_analysis.R, only adding the MCMC-specific fields (iterations, nChains, burnin, thin). The
same pattern works for any other model/data combination already set up for SMC2() in this repo — e.g. EpiSSM_1age on the real
1-age series, or EpiSSM on the simulated multi-age dataset: build opts the way the corresponding SMC2 script already does, add
the MCMC fields, and call PMMH(SSM, opts) instead of SMC2(SSM, opts). 
PMMH() requires the BayesianTools package (install.packages("BayesianTools")); it isn’t needed for any of the SMC2
based scripts.

All three run SMC², plot filtering-distribution and save figures to `figures/real_data/` or `figures/sim_data/`.

**Load order.** `smc2.R` sources `utils.R`, `resampling.R`,
`particle_filter.R`, `priors.R`, and `mcmc_kernels.R` itself, so in most
scripts you only need:

```r
source("src/smc2.R")
source("src/posterior_marginal.R")
source("src/delay_distributions.R")
```

`mcmc_kernels.R` calls `PF_SMC2()` (defined in `smc2.R`), and `smc2.R`
calls the kernels in `mcmc_kernels.R`. This mutual reference is fine in R
— function calls are resolved when they run, not when the file is
sourced — but both files must be sourced before `SMC2()` is called.


## Simulated-data caching

`load_or_simulate()` (in `simulated_data_cache.R`) runs a simulator
once, caches every element of its output to CSV under
`data/simulated/`, and on every later call just reads the CSVs back —
delete the relevant `data/simulated/<label>_*` files to force a fresh
simulation. `run_multi_simulation.R` uses this for its synthetic
dataset.

## Block vs. global resampling

`OneStepUpdate()` / `BootstrapPF()` (and, through it, `IncrementalPF()`
and `SMC2()`) support two resampling scopes via `opts$resample_scope`:

- **`"block"`** (default) — each age group is resampled independently
  from its own observation weights. This is the original behaviour:
  it improves particle efficiency per age group.
- **`"global"`** — a single set of resampling indices, drawn from the
  joint (product-across-age) weight, is applied to the entire state
  vector at once. This preserves full posterior coherence across age
  groups, at the cost of faster particle degeneracy in high dimension.

Both give the same (unbiased) log-likelihood increment; they only
differ in how particle diversity propagates across age groups. Set it
in `opts` before calling `BootstrapPF()` or `SMC2()`:

```r
opts$resample_scope  <- "global"   # or "block" (default)
opts$resample_method <- "systematic"  # "multinomial" | "stratified"
```

