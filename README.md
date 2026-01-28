# latentmed

**Automatic Mediation Analysis Under Measurement Error Via Bayesian Machine Learning**

`latentmed` implements Bayesian nonparametric methods for causal mediation analysis when mediators and/or outcomes are latent constructs measured through multiple indicators. The package uses BART (Bayesian Additive Regression Trees) via [SoftBart](https://github.com/theodds/SoftBart) for flexible, automatic modeling of mediator and outcome regressions.

## Installation

```r
# Install from GitHub
devtools::install_github("XinranSong/latentmed")

# Or install from local source
devtools::install("path/to/latentmed")
```

### Dependencies

```r
install.packages(c("MASS", "SoftBart", "progress"))
```

## Usage

### Basic Example

```r
library(latentmed)

# Prepare data
# M: N x J matrix of mediator indicators (J >= 3 for latent, J = 1 for observed)
# Y: N x K matrix of outcome indicators (K >= 3 for latent, K = 1 for observed)
# X: N x P matrix with treatment in first column

# Latent mediator, observed outcome
results <- run_npgibbs(X, Y, M, treatment_col = 1)

# Both latent mediator and latent outcome
results <- run_npgibbs_latent(X, Y, M, treatment_col = 1)

# Automatic detection (recommended)
results <- run_npgibbs_general(X, Y, M, treatment_col = 1)
```

### Extracting Results

```r
# Direct effects (NDE)
mean(results$direct_0)   # NDE with M fixed at M(0)
mean(results$direct_1)   # NDE with M fixed at M(1)

# Indirect effects (NIE)
mean(results$indirect_0) # NIE with treatment fixed at 0
mean(results$indirect_1) # NIE with treatment fixed at 1

# Total effect
total_effect <- results$direct_0 + results$indirect_0

# 95% Credible intervals
quantile(results$indirect_0, c(0.025, 0.975))

# Population average treatment effects (PATE)
mean(results$indirect_0_pate)
```

## Model

The package implements the following structural model:

**Structural Equations:**
- Mediator: $\theta_M = f_M(X) + \epsilon_M$
- Outcome: $\theta_Y = f_Y(X, \theta_M) + \epsilon_Y$

**Measurement Models:**
- Mediator indicators: $M_j = a_j^M + b_j^M \theta_M + e_j^M$
- Outcome indicators: $Y_k = a_k^Y + b_k^Y \theta_Y + e_k^Y$

Where $f_M$ and $f_Y$ are modeled flexibly using BART.

**Identification:** First indicator loading fixed to 1, intercept to 0.

## Main Functions

| Function | Description |
|----------|-------------|
| `run_npgibbs()` | Gibbs sampler for latent M, observed Y |
| `run_npgibbs_latent()` | Gibbs sampler for latent M, latent Y |
| `run_npgibbs_general()` | Unified sampler (auto-detects latent/observed) |
| `gibbs_factor_analysis()` | Factor analysis for initialization |
| `partial_linear_cma()` | Partial linear model alternative |

## Stan Models

For fully parametric (linear) versions, Stan models are included:

```r
stan_path <- system.file("stan", "latent_cma_yobs.stan", package = "latentmed")
```

## File Structure

```
latentmed/
├── R/
│   ├── np_gibbs_general.R      # Main unified Gibbs sampler
│   ├── npgibbs.R               # Latent M, latent Y
│   ├── npgibbs_yobs.R          # Latent M, observed Y
│   ├── gibbs_factor_analysis.R # Factor analysis sampler
│   ├── partial_linear_cma.R    # Partial linear model
│   ├── em_estimate.R           # EM for initialization
│   ├── set_prior_partial_linear.R
│   ├── nonlinear_simulation_functions.R
│   └── pd_plot2.R              # Partial dependence plots
├── inst/stan/
│   ├── latent_cma.stan
│   └── latent_cma_yobs.stan
└── man/                        # Documentation
```

## License

MIT
