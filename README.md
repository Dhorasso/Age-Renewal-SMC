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
| `utils.R`                 | `%||%`, `logSumExp`, progress bar, `par_lapply` |
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
| `epi_diagnostics.R`       | Next-generation-matrix R_t / R_{a,t}, sort-then-sum population totals |
| `diagnostics.R`           | All ggplot2 plotting helpers (filtering ribbons, SMC² diagnostics, prior/posterior overlays, grid layout) |
| `sim_data_1age.R`         | Generation-time / report-delay kernels for the 1-age driver |
| `sim_data_multi.R`        | Forward-simulates synthetic age-stratified data for the simulation study |
| `simulated_data_cache.R`  | `load_or_simulate()`: cache any simulator's output to CSV under `data/simulated/` |

## `scripts/` (entry points)

| Script | Fits | Data |
|---|---|---|
| `run_1age_covid.R` | 1-age EpiSSM | real total daily cases (`data/real/Covid_Data_Ireland.csv`) |
| `run_multi_simulation.R` | age-stratified EpiSSM | synthetic data, cached to `data/simulated/` |
| `run_ireland_analysis.R` | age-stratified EpiSSM | real age-stratified cases, with 14-day forecast + CRPS |

All three run SMC², plot filtering-distribution estimates only (no PMMH — dropped per project decision), and save figures to `figures/real_data/` or `figures/sim_data/`.

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

## Renaming / reorganisation notes (vs. the original scripts)

- `BPFsmc2.R` → **`particle_filter.R`** (name reflects its content: the
  bootstrap particle filter, used both standalone and inside SMC²).
- `SMC2.R` was a single ~650-line file mixing five concerns. It has been
  split into `smc2.R` (control flow only), `priors.R` (prior
  sampling/density), and `mcmc_kernels.R` (move-step kernels + adaptive
  `Nx`), with shared helpers (`%||%`, `logSumExp`, progress bar,
  `par_lapply`) factored out into `utils.R` instead of being redefined
  in three different files.
- `delay_pmf.R` → **`delay_distributions.R`**; it contained three
  near-duplicate analytic discretizers (`discretize_gamma`,
  `discr_gamma`, `discretize_dist`). These are consolidated into one
  `discretize_dist(dist = "gamma" | "lognormal")`.
- `PosteriorMarginal.R` → **`posterior_marginal.R`** (lowercase, for
  consistent file naming across the package).
- `Resampling.R` → **`resampling.R`** (same reason); the systematic and
  stratified search loop is factored into one shared helper
  (`.cdf_lookup`) instead of being duplicated inline.
- `EpiSSM_1age.R` → **`epi_ssm_1age.R`**; `EpiSSM_multi_data.R` →
  **`epi_ssm_multiage.R`**.
- `Diag_smc2.R` → **`diagnostics.R`**. `plot_filter_one_age()`,
  `.assemble_age_fig()`, and `.grid_shared_ylabel()` were each
  re-defined (with small, inconsistent variations) inline in
  `test_1age_smc2.R`, `Test_muti_data.R`, and `Ireland_analysis.R`.
  They now live only in `diagnostics.R`. All PMMH-boxplot helpers
  (`.add_pmmh_boxplot()`) have been removed.
- The next-generation-matrix / R_t / R_{a,t} / sort-then-sum logic,
  duplicated across `Test_muti_data.R` and `Ireland_analysis.R`, is
  consolidated into **`epi_diagnostics.R`**.
- `test_1age_smc2.R` → **`scripts/run_1age_covid.R`**;
  `Test_muti_data.R` → **`scripts/run_multi_simulation.R`**;
  `Ireland_analysis.R` → **`scripts/run_ireland_analysis.R`**.

## Bugs found and fixed

1. **`epi_ssm_1age.R`, `ObsProcess()`**: when `week_effect = FALSE`, the
   original code set `phi <- theta[2]` but then called
   `rnbinom(..., size = kappa, ...)` — `kappa` was never defined in
   that branch and would error at runtime. Fixed by using `kappa` in
   both branches (confirmed with project owner).
2. **`Diag_smc2.R`, `plot_filter_one_age()`**: the `true_points` branch
   referenced a global `cut_date` that wasn't a parameter of the
   function — it only worked because a global variable of that name
   happened to be defined later in `Ireland_analysis.R`. Replaced with
   a `.anchor_points()` helper that anchors `true_points` at
   `forecast_start_date + 1` when supplied (confirmed with project
   owner), falling back to tail-alignment of the panel's own time axis
   otherwise.
3. **Population-CSV collapsing bug** in the (now-consolidated)
   contact-matrix aggregation code: the original driver scripts always
   collapsed the population vector's upper tail assuming a
   single-year-of-age table running past age 85 (as the original,
   uploaded `PEA11...csv` did). Given `R 86:85` counts *backwards*
   (`c(86, 85)`), running this against
   `Ireland_country_level_age_distribution_85.csv` — which already has
   exactly 85 age bands — would silently corrupt the final population
   entry. `read_age_distribution()` in `contact_matrix_utils.R` now
   only collapses when the file actually has more rows than the
   target resolution.
4. **Inconsistent next-generation-matrix convention**: `Test_muti_data.R`
   computed `K[i,j] = beta_j * C[i,j]` while `Ireland_analysis.R` computed
   `K[i,j] = beta_i * C[i,j]` for the *same* model. Only the latter
   matches the renewal equation actually implemented in
   `epi_ssm_multiage.R::StateProcess()` (age `a`'s new infections are
   driven by its own `beta_a(t)`, applied to contacts received from
   every source age `b`, weighted by `ContMatrix[b, a]`). `epi_diagnostics.R`
   documents this derivation explicitly and uses the correct convention
   throughout — worth double-checking against your own derivation.
5. **`epi_ssm_multiage.R`**: contact-matrix reparameterisation helpers
   (`.contact_row()`, Cholesky/scale/estimate modes, gated by
   `opts$chol_contact` / `scale_contact` / `estimate_contact`) exist but
   are never actually called by `StateProcess()`, which always reads
   `opts$ContMatrix` directly. This looks like unfinished work rather
   than a typo, so it has been **left as-is** rather than guessed at —
   if you intended one of those modes to be active, let me know and
   I'll wire it in.

## Missing dependencies (not part of either upload)

- `src/sim_data1age.R` and `src/sim_multi_data.R` were referenced by
  the original drivers but never uploaded. `sim_data_1age.R` was a
  safe substitute (the 1-age driver only ever used it for two delay
  kernels). `sim_data_multi.R`, however, had to *forward-simulate* an
  internally consistent synthetic dataset from the model itself,
  since the multi-age driver validates against ground-truth
  `beta_true`/`I_true` trajectories that weren't provided — see that
  file's header for the (illustrative, replaceable) parameter choices.
- `src/Marginal_NPF.R` (`get_marginal_states()`/`simulate_observations()`)
  was referenced only by the 1-age driver; the other two drivers
  already used the equivalent `PosteriorMarginal()`, so
  `run_1age_covid.R` has been switched to use that instead.

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

## Package-readiness

Functions use explicit, documented arguments and roxygen-style comment
blocks (`#'`) so that converting `src/` into an installable package
(`devtools::create()` + move files to `R/`, add a `DESCRIPTION` and
`NAMESPACE`) is a small step rather than a rewrite. External
dependencies used so far: `compiler`, `parallel`, `MASS`, `FNN`.
