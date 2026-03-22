
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixtureEM <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**mixtureEM** is an R package for person-centred mixture modelling. It
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
pak::pak("pdvalencia/mixtureEM")
```

## Quick start

### Binary LCA

``` r
library(mixtureEM)

# Simulate binary data from a 3-class population
set.seed(42)
probs <- list(
  c(0.1, 0.2, 0.1, 0.1, 0.1),  # Class 1: Low on all items
  c(0.9, 0.8, 0.7, 0.1, 0.1),  # Class 2: High on items 1–3
  c(0.8, 0.8, 0.7, 0.9, 0.9)   # Class 3: High on all items
)
weights <- c(0.6, 0.3, 0.1)
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
#> Converged          : TRUE (in 11 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -778.16
#>   Rel. Entropy   : 0.8034
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 61.31%
#>   Class 2: 25.02%
#>   Class 3: 13.67%
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
#> Item_1               |   0.138 |   0.873 |   0.783
#> Item_2               |   0.248 |   0.894 |   0.834
#> Item_3               |   0.116 |   0.774 |   0.893
#> Item_4               |   0.131 |   0.027 |   0.851
#> Item_5               |   0.071 |   0.081 |   0.873
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
#> Assigned Class 1    0.966    0.027    0.007
#> Assigned Class 2    0.130    0.817    0.053
#> Assigned Class 3    0.036    0.052    0.912
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
#>   Classes       LL Params      AIC      BIC    SABIC Entropy
#> 1       1 -905.532      5 1821.063 1839.582 1823.725   1.000
#> 2       2 -805.129     11 1632.257 1672.999 1638.113   0.767
#> 3       3 -778.180     17 1590.360 1653.324 1599.410   0.807
#> 4       4 -777.025     23 1600.050 1685.237 1612.294   0.662
#> 
#> -> Best model according to BIC: 3 classes

# The best-fitting model and its full table
sel$best_k
#> [1] 3
sel$fit_table
#>   Classes        LL Params      AIC      BIC    SABIC   Entropy
#> 1       1 -905.5316      5 1821.063 1839.582 1823.725 1.0000000
#> 2       2 -805.1285     11 1632.257 1672.999 1638.113 0.7671238
#> 3       3 -778.1798     17 1590.360 1653.324 1599.410 0.8073438
#> 4       4 -777.0249     23 1600.050 1685.237 1612.294 0.6624666
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
#>   Intercept         0.779  [ 0.544,  1.116]     0.173
#>   z1                1.976  [ 1.483,  2.633]     0.000
#> 
#> Class 2 ON
#>   Intercept         0.906  [ 0.654,  1.256]     0.555
#>   z1                0.417  [ 0.307,  0.565]     0.000
#> =========================================================
```

``` r
# Odds ratios relative to Class 3 (the reference)
coef(fit_cov, ref_class = 3)
#>               Intercept        z1
#> Class 1       0.7791112 1.9761065
#> Class 2       0.9062827 0.4166476
#> Class 3 (Ref) 1.0000000 1.0000000
```

``` r
# Confidence intervals for odds ratios
confint(fit_cov, ref_class = 3)
#> Calculating Analytical Confidence Intervals (Hessian)...
#> $Intercept
#>                  OR Lower Upper
#> Class 1       0.779 0.544 1.116
#> Class 2       0.906 0.654 1.256
#> Class 3 (Ref) 1.000 1.000 1.000
#> 
#> $z1
#>                  OR Lower Upper
#> Class 1       1.976 1.483 2.633
#> Class 2       0.417 0.307 0.565
#> Class 3 (Ref) 1.000 1.000 1.000
```

### Hypothesis tests for covariates

#### Analytical Wald test

Fast, no resampling required — uses the observed information matrix.

``` r
analytical_wald_test(fit_cov, term_name = "z1", ref_class = 3)
#>   Covariate Wald_Chi2 df p_value     Method
#> 1        z1     73.45  2       0 Analytical
```

#### Bootstrap Wald omnibus test

More robust in small samples or when classes overlap.

``` r
boot <- bootstrap_covariates(fit_cov, X, Z, n_reps = 200, ref_class = 3)
#> Running bootstrap with alignment (Ref Class: 3, 200 reps)...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
boot$p_values
#>               Intercept        z1
#> Class 1       0.2525765 0.1286137
#> Class 2       0.6574223 0.1476802
#> Class 3 (Ref)       NaN       NaN
wald_omnibus_test(boot, term_name = "z1")
#>   Covariate Wald_Chi2 df p_value
#> 1        z1     4.405  2  0.1105
confint(fit_cov, boot_results = boot, ref_class = 3)
#> Calculating Bootstrapped Confidence Intervals...
#> $Intercept
#>                  OR Lower Upper
#> Class 1       0.779 0.508 1.195
#> Class 2       0.906 0.587 1.400
#> Class 3 (Ref) 1.000 1.000 1.000
#> 
#> $z1
#>                  OR Lower Upper
#> Class 1       1.976 0.821 4.757
#> Class 2       0.417 0.127 1.363
#> Class 3 (Ref) 1.000 1.000 1.000
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
#> Running BLRT: 2 vs 3 classes (This may take a while)...
#>   Bootstrap draw 10 / 100 completed.
#>   Bootstrap draw 20 / 100 completed.
#>   Bootstrap draw 30 / 100 completed.
#>   Bootstrap draw 40 / 100 completed.
#>   Bootstrap draw 50 / 100 completed.
#>   Bootstrap draw 60 / 100 completed.
#>   Bootstrap draw 70 / 100 completed.
#>   Bootstrap draw 80 / 100 completed.
#>   Bootstrap draw 90 / 100 completed.
#>   Bootstrap draw 100 / 100 completed.
result$p_value
#> [1] 0.00990099
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

If you use mixtureEM in published research, please cite it as:

    Valencia, P. D. (2025). mixtureEM: Latent Class and Profile Analysis in R.
    R package. https://github.com/pdvalencia/mixtureEM
