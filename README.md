
<!-- README.md is generated from README.Rmd. Please edit that file -->

# unmixR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/unmixR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/unmixR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**unmixR** is an R package for person-centred mixture modelling. It
provides a unified, flexible interface for **latent class analysis
(LCA)** and **latent profile analysis (LPA)**, supporting a wide range
of measurement models, structural models, and multi-step estimation
strategies.

## Features

- **Measurement models:** binary (Bernoulli), categorical (Multinoulli),
  and continuous (Gaussian) indicators, including missing-data variants
- **Mixed measurement models:** combine different variable types in a
  single model via a named-list descriptor
- **Structural models:** covariate predictors of class membership, and
  distal outcomes (continuous, categorical — pooled or class-specific
  regression)
- **Multi-step estimation:** 1-step (simultaneous), 2-step, and 3-step
  with **BCH** or **ML** bias correction
- **Model selection:** AIC, BIC, SABIC, and relative entropy across a
  range of K
- **Bootstrap Likelihood Ratio Test (BLRT)** for comparing nested models
- **Inference for covariates:** analytical Wald tests, bootstrap
  standard errors, and confidence intervals for odds ratios
- **Survey weights** support throughout

## Installation

``` r
# install.packages("pak")
pak::pak("pdvalencia/unmixR")
```

## Quick start

### Binary LCA

``` r
library(unmixR)

# Simulate binary data from a 3-class population
set.seed(42)
probs <- list(
  c(0.9, 0.8, 0.7, 0.1, 0.1),  # Class 1: high on items 1–3
  c(0.1, 0.2, 0.1, 0.9, 0.8),  # Class 2: high on items 4–5
  c(0.5, 0.5, 0.5, 0.5, 0.5)   # Class 3: uniform
)
weights <- c(0.4, 0.4, 0.2)
n <- 300
classes <- sample(1:3, n, replace = TRUE, prob = weights)
X <- t(sapply(classes, function(k) rbinom(5, 1, probs[[k]])))

# Fit a 3-class LCA
fit <- fit_mixture(X,
                   n_components  = 3,
                   measurement   = "binary",
                   n_init        = 5,
                   random_state  = 42)
print(fit)
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 3
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 6 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -876.49
#>   Rel. Entropy   : 0.5452
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 43.87%
#>   Class 2: 30.97%
#>   Class 3: 25.17%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```

``` r
# Item-response probabilities for each class
measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2 | Class 3
#> -------------------------------------------------- 
#> Item_1               |   0.892 |   0.151 |   0.232
#> Item_2               |   0.799 |   0.288 |   0.292
#> Item_3               |   0.746 |   0.155 |   0.245
#> Item_4               |   0.076 |   0.739 |   0.759
#> Item_5               |   0.083 |   0.779 |   0.752
#> =========================================================
```

``` r
# Average posterior probabilities — diagonal values near 1 = well-separated classes
classification_diagnostics(fit)
#> =========================================================
#>           AVERAGE POSTERIOR PROBABILITIES (AvePP)        
#> =========================================================
#> Rows: Modal Assignment | Columns: Mean Probability
#> 
#>                  Prob C 1 Prob C 2 Prob C 3
#> Assigned Class 1    0.925    0.034    0.041
#> Assigned Class 2    0.012    0.594    0.393
#> Assigned Class 3    0.057    0.427    0.516
#> =========================================================
```

### Model selection

Use `compare_mixtures()` to fit K = 1 to 4 and compare fit indices. BIC
is the most common criterion for determining the number of classes.

``` r
sel <- compare_mixtures(X,
                        k_range     = 1:4,
                        measurement = "binary",
                        n_init      = 5)
#> Running Model Selection across K = 1 to 4...
#> 
#> Fitting 1-class model...
#> Fitting 2-class model...
#> Fitting 3-class model...
#> Fitting 4-class model...
#> 
#> === Model Selection Summary ===
#>   Classes        LL Params      AIC      BIC    SABIC Entropy
#> 1       1 -1035.219      5 2080.437 2098.956 2083.099   1.000
#> 2       2  -876.714     11 1775.428 1816.169 1781.284   0.820
#> 3       3  -876.118     17 1786.235 1849.199 1795.286   0.573
#> 4       4  -876.071     23 1798.141 1883.328 1810.386   0.455
#> 
#> -> Best model according to BIC: 2 classes

# The best-fitting model and its full table
sel$best_k
#> [1] 2
sel$fit_table
#>   Classes         LL Params      AIC      BIC    SABIC   Entropy
#> 1       1 -1035.2185      5 2080.437 2098.956 2083.099 1.0000000
#> 2       2  -876.7139     11 1775.428 1816.169 1781.284 0.8195552
#> 3       3  -876.1176     17 1786.235 1849.199 1795.286 0.5733766
#> 4       4  -876.0706     23 1798.141 1883.328 1810.386 0.4550529
```

### Latent Profile Analysis (LPA)

For continuous indicators, switch `measurement` to `"continuous"`:

``` r
# Simulate continuous data from a 2-profile population
set.seed(1)
n1 <- 150; n2 <- 150
X_cont <- rbind(
  matrix(rnorm(n1 * 4, mean = c(2, 2, -2, -2), sd = 1), nrow = n1, byrow = TRUE),
  matrix(rnorm(n2 * 4, mean = c(-2, -2, 2, 2), sd = 1), nrow = n2, byrow = TRUE)
)

fit_lpa <- fit_mixture(X_cont,
                       n_components = 2,
                       measurement  = "continuous",
                       n_init       = 5,
                       random_state = 1)
print(fit_lpa)
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 6 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -1939.80
#>   Rel. Entropy   : 1.0000
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 50.00%
#>   Class 2: 50.00%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
measurement_summary(fit_lpa)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CONTINUOUS MEANS
#> Indicator            | Class 1 | Class 2
#> ---------------------------------------- 
#> Item_1               |   1.949 |  -2.057
#> Item_2               |   2.085 |  -2.101
#> Item_3               |  -1.996 |   2.012
#> Item_4               |  -1.991 |   1.915
#> =========================================================
```

### 3-step LCA with a covariate

The 3-step approach first estimates the classes on the measurement model
alone, then regresses class membership on external covariates — avoiding
the parameter contamination that occurs in 1-step estimation. The `"ML"`
correction adjusts for classification error.

``` r
# Covariate correlated with class membership
set.seed(7)
Z <- matrix(ifelse(classes == 1, rnorm(n, 1, 1),
            ifelse(classes == 2, rnorm(n, -1, 1),
                                  rnorm(n,  0, 1))),
            ncol = 1)
colnames(Z) <- "z1"

fit_cov <- fit_mixture(X, Y = Z,
                       n_components = 3,
                       measurement  = "binary",
                       structural   = "covariate",
                       n_steps      = 3,
                       correction   = "ML",
                       n_init       = 5,
                       random_state = 7)

# Structural model: odds ratios with 95 % CIs (analytical SEs)
summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 3
#> ---------------------------------------------------------
#>                      OR       [95% CI]        P-Value
#> 
#> Class 1 ON
#>   Intercept         1.110  [ 0.786,  1.569]     0.553
#>   z1                0.548  [ 0.406,  0.741]     0.000
#> 
#> Class 2 ON
#>   Intercept         0.644  [ 0.414,  1.002]     0.051
#>   z1                5.325  [ 3.361,  8.439]     0.000
#> =========================================================
```

``` r
# Odds ratios relative to Class 3 (the reference)
coef(fit_cov, ref_class = 3)
#>               Intercept       z1
#> Class 1       1.1104436 0.548481
#> Class 2       0.6439104 5.325464
#> Class 3 (Ref) 1.0000000 1.000000
```

``` r
# Confidence intervals for odds ratios
confint(fit_cov, ref_class = 3)
#> Calculating Analytical Confidence Intervals (Hessian)...
#> $Intercept
#>                  OR Lower Upper
#> Class 1       1.110 0.786 1.569
#> Class 2       0.644 0.414 1.002
#> Class 3 (Ref) 1.000 1.000 1.000
#> 
#> $z1
#>                  OR Lower Upper
#> Class 1       0.548 0.406 0.741
#> Class 2       5.325 3.361 8.439
#> Class 3 (Ref) 1.000 1.000 1.000
```

### Hypothesis tests for covariates

#### Analytical Wald test

Fast, no resampling required — uses the observed information matrix.

``` r
analytical_wald_test(fit_cov, term_name = "z1", ref_class = 3)
#>   Covariate Wald_Chi2 df p_value     Method
#> 1        z1    81.155  2       0 Analytical
```

#### Bootstrap Wald omnibus test

More robust in small samples or when classes overlap.

``` r
boot <- bootstrap_covariates(fit_cov, X, Z, n_reps = 200, ref_class = 3)
boot$p_values
wald_omnibus_test(boot, term_name = "z1")
confint(fit_cov, boot_results = boot, ref_class = 3)
```

### Bootstrap Likelihood Ratio Test (BLRT)

Tests whether a K-class model fits significantly better than a
(K−1)-class model, using a parametric bootstrap to approximate the null
distribution.

``` r
result <- calc_blrt(X,
                    k_small      = 2,
                    k_large      = 3,
                    measurement  = "binary",
                    n_reps       = 100)
result$p_value
```

## Function reference

| Function                       | Description                                                    |
|--------------------------------|----------------------------------------------------------------|
| `fit_mixture()`                | Fit an LCA / LPA model (core function)                         |
| `compare_mixtures()`           | Fit models across a range of K and compare fit indices         |
| `print()`                      | Print a compact model summary                                  |
| `summary()`                    | Print structural model parameters (OR, CI, p-values)           |
| `measurement_summary()`        | Print item-response probabilities or means per class           |
| `classification_diagnostics()` | Average posterior probability (AvePP) matrix                   |
| `coef()`                       | Extract odds ratios from a covariate model                     |
| `confint()`                    | Confidence intervals for odds ratios (analytical or bootstrap) |
| `analytical_wald_test()`       | Analytical Wald chi-squared test for a covariate               |
| `bootstrap_covariates()`       | Bootstrap standard errors for covariate parameters             |
| `wald_omnibus_test()`          | Bootstrap Wald omnibus test for a covariate                    |
| `calc_blrt()`                  | Bootstrap Likelihood Ratio Test for class enumeration          |

## Supported model types

### Measurement (`measurement =`)

| String                                     | Indicator type                         |
|--------------------------------------------|----------------------------------------|
| `"binary"` / `"bernoulli"`                 | Binary (0/1)                           |
| `"binary_nan"` / `"bernoulli_nan"`         | Binary with missing data               |
| `"categorical"` / `"multinoulli"`          | Ordinal / polytomous                   |
| `"categorical_nan"` / `"multinoulli_nan"`  | Ordinal with missing data              |
| `"continuous"` / `"gaussian_diag"`         | Continuous (estimated variance)        |
| `"continuous_nan"` / `"gaussian_diag_nan"` | Continuous with missing data           |
| `"gaussian"` / `"gaussian_unit"`           | Continuous (unit variance)             |
| Named list                                 | Mixed model combining any of the above |

### Structural (`structural =`)

| String                           | Model                                                 |
|----------------------------------|-------------------------------------------------------|
| `"covariate"`                    | Multinomial logistic regression predicting class      |
| `"distal_continuous"`            | Continuous distal outcome (class means)               |
| `"distal_continuous_regression"` | Continuous distal outcome moderated by a covariate    |
| `"distal_pooled"`                | Categorical distal outcome, pooled regression         |
| `"distal_regression"`            | Categorical distal outcome, class-specific regression |

## Citation

If you use unmixR in published research, please cite it as:

    Valencia, P. (2025). unmixR: Latent Class and Profile Analysis in R.
    R package. https://github.com/pdvalencia/unmixR
