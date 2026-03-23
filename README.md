
# mixtureEM <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml)
**mixtureEM** is an R package for person-centred mixture modelling. It
provides a unified, flexible interface for **latent class analysis
(LCA)** and **latent profile analysis (LPA)**, supporting a wide range
of measurement models, structural models, and multi-step estimation
strategies.

## Features

- **Measurement models:** binary (Bernoulli), categorical (Multinoulli),
  and continuous (Gaussian) indicators, including missing-data variants.
- **Mixed measurement models:** combine different variable types in a
  single model via a named-list descriptor.
- **Structural models:** covariate predictors of class membership, and
  distal outcomes (continuous or categorical — with options for pooled
  or class-specific regression slopes).
- **Multi-step estimation:** 1-step (simultaneous), 2-step, and 3-step
  with **BCH** or **ML** bias correction.
- **Model selection:** AIC, BIC, SABIC, and relative entropy across a
  range of K.
- **Bootstrap Likelihood Ratio Test (BLRT)** for comparing nested
  models.
- **Inference for covariates:** analytical Wald tests, bootstrap
  standard errors, and confidence intervals for odds ratios.
- **Survey weights** support throughout.

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
n1 <- 180; n2 <- 120
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
#> Converged          : TRUE (in 3 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -1935.52
#>   Rel. Entropy   : 1.0000
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 59.97%
#>   Class 2: 40.03%
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
#> Item_1               |   1.920 |  -2.014
#> Item_2               |   2.036 |  -2.076
#> Item_3               |  -1.997 |   2.015
#> Item_4               |  -2.025 |   1.942
#> =========================================================
```

### 3-step LCA with a covariate

The 3-step approach first estimates the classes on the measurement model
alone, then regresses class membership on external covariates — avoiding
the parameter contamination that occurs in 1-step estimation. The `"ML"`
correction is recommended for covariate structural models (Vermunt,
2010).

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

# Structural model: odds ratios with 95% CIs (analytical SEs)
summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> ---------------------------------------------------------
#>                      OR       [95% CI]        P-Value
#> 
#> Class 2 ON
#>   Intercept         1.163  [ 0.802,  1.688]     0.426
#>   z1                0.211  [ 0.148,  0.301]     0.000
#> 
#> Class 3 ON
#>   Intercept         1.284  [ 0.896,  1.838]     0.173
#>   z1                0.506  [ 0.380,  0.674]     0.000
#> =========================================================

# Changing the reference class to 3
summary(fit_cov, ref_class = 3)
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

### Hypothesis tests for covariates

#### Analytical Wald test

Fast, no resampling required — uses the observed information matrix.

``` r
analytical_wald_test(fit_cov, term_name = "z1", ref_class = 1)
#>   Covariate Wald_Chi2 df p_value     Method
#> 1        z1     73.45  2       0 Analytical
```

#### Bootstrap Wald omnibus test

Experimental. Might be more robust in small samples or when classes
overlap.

``` r
boot <- bootstrap_covariates(fit_cov, X, Z, n_reps = 200, ref_class = 1)
#> Running bootstrap with alignment (Ref Class: 1, 200 reps)...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
boot$p_values
#>               Intercept           z1
#> Class 1 (Ref)       NaN          NaN
#> Class 2       0.5952306 4.713721e-07
#> Class 3       0.2525765 1.286137e-01
wald_omnibus_test(boot, term_name = "z1")
#>   Covariate Wald_Chi2 df p_value
#> 1        z1    27.687  2       0
confint(fit_cov, boot_results = boot, ref_class = 3)
#> Calculating Bootstrapped Confidence Intervals...
#> $Intercept
#>                  OR Lower Upper
#> Class 1       0.779 0.779 0.779
#> Class 2       0.906 0.519 1.583
#> Class 3 (Ref) 1.000 0.652 1.534
#> 
#> $z1
#>                  OR Lower Upper
#> Class 1       1.976 1.976 1.976
#> Class 2       0.417 0.227 0.763
#> Class 3 (Ref) 1.000 0.415 2.407
```

### Distal outcomes

A **distal outcome** is a variable that is caused by (rather than
constitutive of) the latent classes. The 3-step approach first
establishes the class solution from the measurement model alone, then
links it to the outcome — protecting the class solution from
contamination.

The choice of bias correction depends on the outcome type:

- **BCH** is recommended for continuous distal outcomes
  (`distal_continuous`, `distal_continuous_pooled`,
  `distal_continuous_regression`). It re-weights observations by the
  inverse of the classification-error matrix without requiring
  distributional assumptions about the outcome (Bakk & Vermunt, 2016).
- **ML** is recommended for covariate and categorical distal models. It
  iterates jointly over measurement and structural components and
  handles missing data on the outcome (Vermunt, 2010).

All distal models share the same `Y` column convention: **column 1 is
always the outcome variable**, and any additional columns are treated as
covariates.

#### Continuous distal: Outcome ~ Class

Use `structural = "distal_continuous"` to estimate the mean of a
continuous outcome in each latent class, with no covariate adjustment.
`summary()` reports class-specific means with 95% CIs and robust
sandwich SEs.

``` r
# Simulate a continuous distal outcome with class-specific means
set.seed(21)
D <- matrix(ifelse(classes == 1, rnorm(n,  0.0, 1),
            ifelse(classes == 2, rnorm(n,  2.5, 1),
                                 rnorm(n, -1.5, 1))),
            ncol = 1)
colnames(D) <- "outcome"

fit_distal <- fit_mixture(X, Y = D,
                          n_components = 3,
                          measurement  = "binary",
                          structural   = "distal_continuous",
                          n_steps      = 3,
                          correction   = "BCH",
                          n_init       = 5,
                          random_state = 21)

# Class-specific means with 95% CIs and robust SEs
summary(fit_distal)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CONTINUOUS DISTAL OUTCOME (MEANS)
#> ---------------------------------------------------------
#>                  Mean       [95% CI]        SE
#>   Class 1        0.050  [-0.096,  0.196]     0.074
#>   Class 2        2.409  [ 2.110,  2.707]     0.152
#>   Class 3       -0.231  [-0.902,  0.441]     0.343
#> =========================================================
```

#### Continuous distal with covariate: Outcome ~ Class + Covariate, pooled slope

Use `structural = "distal_continuous_pooled"` to estimate class-varying
intercepts but a **single covariate slope shared across all classes**.
This is a parsimonious model assuming the covariate’s effect on the
outcome is the same regardless of class membership. The `Y` matrix must
have the outcome in column 1 and covariate(s) in columns 2+.

``` r
# Y layout: outcome (col 1), covariate (col 2)
Y_cont_reg <- cbind(D, Z)
colnames(Y_cont_reg) <- c("outcome", "z1")

fit_distal_pool_cont <- fit_mixture(X, Y = Y_cont_reg,
                                    n_components = 3,
                                    measurement  = "binary",
                                    structural   = "distal_continuous_pooled",
                                    n_steps      = 3,
                                    correction   = "BCH",
                                    n_init       = 5,
                                    random_state = 21)

# Class-specific intercepts + one pooled covariate slope
summary(fit_distal_pool_cont)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CONTINUOUS DISTAL POOLED REGRESSION (Main Effects)
#> ---------------------------------------------------------
#> 
#>   Latent Class (Intercepts):
#>                  Estimate   [95% CI]        P-Value
#>     Class 1        0.294  [ 0.081,  0.506]     0.007
#>     Class 2        2.217  [ 1.889,  2.546]     0.000
#>     Class 3       -0.316  [-0.940,  0.308]     0.321
#> 
#>   Covariates (Pooled Slopes):
#>     Z1           -0.271  [-0.435, -0.108]     0.001
#> =========================================================
```

#### Continuous distal with covariate: Outcome ~ Class × Covariate, class-specific slopes

Use `structural = "distal_continuous_regression"` to regress a
continuous outcome on class membership and a covariate simultaneously,
with **class-specific intercepts and slopes** — that is, the covariate
effect is allowed to differ across classes.

``` r
fit_distal_reg <- fit_mixture(X, Y = Y_cont_reg,
                              n_components = 3,
                              measurement  = "binary",
                              structural   = "distal_continuous_regression",
                              n_steps      = 3,
                              correction   = "BCH",
                              n_init       = 5,
                              random_state = 21)

# Class-specific intercepts and slopes (raw coefficients)
summary(fit_distal_reg)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CONTINUOUS DISTAL REGRESSION (Y ~ Z * Class)
#> ---------------------------------------------------------
#> 
#> Class 1:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept       0.112  [-0.098,  0.323]     0.295
#>   Z1             -0.069  [-0.224,  0.086]     0.380
#> 
#> Class 2:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept       2.238  [ 1.821,  2.654]     0.000
#>   Z1             -0.243  [-0.573,  0.088]     0.150
#> 
#> Class 3:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept      -0.550  [-1.121,  0.020]     0.059
#>   Z1             -1.019  [-1.480, -0.559]     0.000
#> =========================================================
```

#### Categorical distal: Outcome ~ Class + Covariate, pooled slope

For a categorical (ordinal or nominal) distal outcome, mixtureEM
provides the same flexibility regarding covariate slopes.

`structural = "distal_pooled"` estimates class-varying intercepts but a
**single covariate slope shared across all classes**. The outcome is
passed in column 1 and the covariate in column 2.

``` r
# Simulate a binary categorical distal outcome (0/1 → internally 1-indexed)
set.seed(33)
log_odds <- ifelse(classes == 1, -1.0,
            ifelse(classes == 2,  1.5, 0.0)) + 0.4 * Z[, 1]
D_cat <- matrix(rbinom(n, 1, plogis(log_odds)), ncol = 1)
colnames(D_cat) <- "cat_outcome"

# Y layout: categorical outcome (col 1), covariate (col 2)
Y_pooled <- cbind(D_cat, Z)
colnames(Y_pooled) <- c("cat_outcome", "z1")

fit_pooled <- fit_mixture(X, Y = Y_pooled,
                          n_components = 3,
                          measurement  = "binary",
                          structural   = "distal_pooled",
                          n_steps      = 3,
                          correction   = "ML",
                          n_init       = 5,
                          random_state = 33)

# Class-specific intercepts + one pooled covariate slope
summary(fit_pooled)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> POOLED DISTAL REGRESSION (MAIN EFFECTS)
#> Reference Class: 1
#> ---------------------------------------------------------
#> 
#> Outcome Category 2 (vs Category 1) ON
#>                      OR       [95% CI]        P-Value
#> 
#>   Latent Class (Main Effect):
#>     Class 2          8.983  [ 4.254, 18.967]     0.000
#>     Class 3          4.533  [ 2.087,  9.848]     0.000
#> 
#>   Covariates (Main Effect):
#>     Z1              1.551  [ 1.233,  1.951]     0.000
#> =========================================================
```

The summary reports odds ratios for the class main effects (contrasted
against `ref_class`) and a single odds ratio for the covariate that
applies equally across all classes.

#### Categorical distal: Outcome ~ Class × Covariate, class-specific slopes

`structural = "distal_regression"` relaxes the pooling assumption by
estimating **class-specific intercepts and slopes** — the covariate
effect on the outcome is free to vary across classes. Use this model
when you expect the covariate to be a moderator of the class–outcome
relationship.

``` r
# Same data as above; now allow the covariate slope to vary by class
fit_moderated <- fit_mixture(X, Y = Y_pooled,
                             n_components = 3,
                             measurement  = "binary",
                             structural   = "distal_regression",
                             n_steps      = 3,
                             correction   = "ML",
                             n_init       = 5,
                             random_state = 33)

# Separate intercept and covariate slope for each class
summary(fit_moderated)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> SIMULTANEOUS DISTAL REGRESSION (MODERATED BY CLASS)
#> ---------------------------------------------------------
#> 
#> Class 1:
#>                      OR       [95% CI]        P-Value
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       0.353  [ 0.191,  0.651]     0.001
#>     Z1              1.554  [ 0.935,  2.583]     0.089
#> 
#> Class 2:
#>                      OR       [95% CI]        P-Value
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       3.564  [ 2.266,  5.603]     0.000
#>     Z1              1.772  [ 1.293,  2.428]     0.000
#> 
#> Class 3:
#>                      OR       [95% CI]        P-Value
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       1.505  [ 0.795,  2.849]     0.209
#>     Z1              1.307  [ 0.827,  2.067]     0.252
#> =========================================================
```

Comparing `fit_pooled` and `fit_moderated` on AIC / BIC (accessible via
`fit$metrics`) helps determine whether allowing class-specific slopes
meaningfully improves fit:

``` r
cat("Pooled:    BIC =", round(fit_pooled$metrics$bic, 2), "\n")
#> Pooled:    BIC = 2061.59
cat("Moderated: BIC =", round(fit_moderated$metrics$bic, 2), "\n")
#> Moderated: BIC = 2072.18
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

| Function                       | Description                                                                                     |
|--------------------------------|-------------------------------------------------------------------------------------------------|
| `fit_mixture()`                | Fit an LCA / LPA model (core function)                                                          |
| `compare_mixtures()`           | Fit models across a range of K and compare fit indices                                          |
| `print()`                      | Print a compact model summary                                                                   |
| `summary()`                    | Print structural model parameters (OR, CI, p-values; or distal means / regression coefficients) |
| `measurement_summary()`        | Print item-response probabilities or means per class                                            |
| `classification_diagnostics()` | Average posterior probability (AvePP) matrix                                                    |
| `coef()`                       | Extract odds ratios from a covariate model                                                      |
| `confint()`                    | Confidence intervals for odds ratios (analytical or bootstrap)                                  |
| `analytical_wald_test()`       | Analytical Wald chi-squared test for a covariate                                                |
| `bootstrap_covariates()`       | Bootstrap standard errors for covariate parameters                                              |
| `wald_omnibus_test()`          | Bootstrap Wald omnibus test for a covariate                                                     |
| `calc_blrt()`                  | Bootstrap Likelihood Ratio Test for class enumeration                                           |

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

| String                           | Outcome type       | Covariate slope           | `Y` layout                            | Recommended correction |
|----------------------------------|--------------------|---------------------------|---------------------------------------|------------------------|
| `"covariate"`                    | — (predicts class) | —                         | Covariates only                       | `"ML"`                 |
| `"distal_continuous"`            | Continuous         | None                      | Outcome only                          | `"BCH"`                |
| `"distal_continuous_pooled"`     | Continuous         | **Pooled** across classes | Outcome (col 1), covariates (cols 2+) | `"BCH"`                |
| `"distal_continuous_regression"` | Continuous         | Class-specific            | Outcome (col 1), covariates (cols 2+) | `"BCH"`                |
| `"distal_pooled"`                | Categorical        | **Pooled** across classes | Outcome (col 1), covariates (cols 2+) | `"ML"`                 |
| `"distal_regression"`            | Categorical        | Class-specific            | Outcome (col 1), covariates (cols 2+) | `"ML"`                 |

> **On the pooled vs. moderated distinction:** mixtureEM provides both
> pooled (main effects only) and moderated (class-specific slopes)
> options for **both categorical and continuous** distal outcomes.
> Choose the pooled version for a more parsimonious model, or the
> regression/moderated version if you suspect the covariate’s effect
> differs across classes.

## Acknowledgements

**mixtureEM** draws strong inspiration from the Python package
[**StepMix**](https://github.com/Labo-Lacourse/stepmix) (Morin et al.,
2025), which pioneered open-source, bias-adjusted multi-step estimation
of generalised mixture models with external variables. The design of the
stepwise estimators, BCH and ML corrections, and the overall modular
measurement–structural model architecture in this package follow the
framework laid out in StepMix. If you make use of these methods, please
also consider citing the StepMix paper:

> Morin, S., Legault, R., Laliberté, F., Bakk, Z., Giguère, C.-É., de la
> Sablonnière, R., & Lacourse, É. (2025). StepMix: A Python package for
> pseudo-likelihood estimation of generalised mixture models with
> external variables. *Journal of Statistical Software*, *113*(8), 1–39.
> <https://doi.org/10.18637/jss.v113.i08>

## Citation

If you use mixtureEM in published research, please cite it as:

    Valencia, P. D. (2026). mixtureEM: Latent Class and Profile Analysis in R.
    R package. [https://github.com/pdvalencia/mixtureEM](https://github.com/pdvalencia/mixtureEM)
