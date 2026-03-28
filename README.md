
# mixtureEM <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml)

**mixtureEM** is an R package for person-centred mixture modelling. It
provides a unified, flexible interface for **latent class analysis
(LCA)** and **latent profile analysis (LPA)**, supporting a wide range
of measurement models, structural models, and multi-step estimation
strategies.

The core idea behind these models is that an unobserved (latent)
categorical variable — the *class* or *profile* — explains patterns of
similarity among observed indicators. LCA is used when indicators are
binary or categorical; LPA when they are continuous.

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

The most common use case: binary (0/1) indicator items, such as symptom
checklists or yes/no questionnaire responses.

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
n       <- 300
classes <- sample(1:3, n, replace = TRUE, prob = weights)
X       <- t(sapply(classes, function(k) rbinom(5, 1, probs[[k]])))

# Fit a 3-class LCA
fit <- fit_mixture(X,
                   n_components = 3,
                   measurement  = "binary",
                   n_init       = 5,
                   random_state = 42)
print(fit)
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 3
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 11 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -777.57
#>   Rel. Entropy   : 0.8236
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 62.28%
#>   Class 2: 24.67%
#>   Class 3: 13.05%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```

``` r
# Item-response probabilities for each class
# Each value is the probability of endorsing that item given class membership
measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2 | Class 3
#> -------------------------------------------------- 
#> Item_1               |   0.140 |   0.888 |   0.793
#> Item_2               |   0.253 |   0.899 |   0.850
#> Item_3               |   0.123 |   0.774 |   0.915
#> Item_4               |   0.133 |   0.006 |   0.912
#> Item_5               |   0.076 |   0.104 |   0.847
#> =========================================================
```

``` r
# Average posterior probabilities — diagonal values near 1 indicate
# well-separated, clearly-defined classes
classification_diagnostics(fit)
#> =========================================================
#>           AVERAGE POSTERIOR PROBABILITIES (AvePP)        
#> =========================================================
#> Rows: Modal Assignment | Columns: Mean Probability
#> 
#>                  Prob C 1 Prob C 2 Prob C 3
#> Assigned Class 1    0.971    0.024    0.005
#> Assigned Class 2    0.141    0.824    0.034
#> Assigned Class 3    0.039    0.010    0.951
#> =========================================================
```

### Model selection

In practice, the number of classes is unknown and must be determined
from the data. `compare_mixtures()` fits models from K = 1 to K = 4 and
returns standard fit indices. **BIC** is the most widely used criterion:
lower is better, and you are looking for the model where BIC stops
meaningfully decreasing.

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
#> 2       2 -805.043     11 1632.086 1672.828 1637.942   0.769
#> 3       3 -777.566     17 1589.133 1652.097 1598.183   0.824
#> 4       4 -775.649     23 1597.298 1682.485 1609.542   0.764
#> 
#> -> Best model according to BIC: 3 classes

# sel$best_k  — the K with the lowest BIC, as an integer for use in subsequent code
# sel$models  — the full list of fitted model objects, indexed as "K1", "K2", etc.
# sel$fit_table — data frame with LL, AIC, BIC, SABIC, and Entropy for each K
```

### Latent Profile Analysis (LPA)

For continuous indicators (e.g., scale scores, physiological measures),
switch `measurement` to `"continuous"`. The model estimates a mean and
variance for each indicator within each profile.

``` r
# Simulate continuous data from a 2-profile population
set.seed(1)
n1 <- 180; n2 <- 120
X_cont <- rbind(
  matrix(rnorm(n1 * 4, mean = c( 2,  2, -2, -2), sd = 1), nrow = n1, byrow = TRUE),
  matrix(rnorm(n2 * 4, mean = c(-2, -2,  2,  2), sd = 1), nrow = n2, byrow = TRUE)
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

**Why 3-step?** The 3-step approach separates the measurement and
structural parts of the model. In the 1-step approach, the covariate
influences how the classes are formed, which conflates class structure
with class predictors. The 3-step approach first establishes the class
solution from the indicators alone (step 1), assigns each person to a
class (step 2), then regresses class membership on the covariate (step
3) — while correcting for the classification error introduced in step 2.
This protects the substantive meaning of the classes.

The `"ML"` correction is recommended for covariate structural models
(Vermunt, 2010).

``` r
# Covariate correlated with class membership
set.seed(7)
Z <- matrix(ifelse(classes == 1, rnorm(n,  1, 1),
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

# Structural model: odds ratios with 95% CIs and p-values
# The reference class (here Class 1) has OR = 1 by definition
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
#>   Intercept         0.371  [ 0.251,  0.549]    < .001
#>   z1                0.158  [ 0.104,  0.241]    < .001
#> 
#> Class 3 ON
#>   Intercept         0.308  [ 0.210,  0.452]    < .001
#>   z1                0.334  [ 0.227,  0.492]    < .001
#> =========================================================

# Change the reference class to Class 3
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
#>   Intercept         3.248  [ 2.214,  4.765]    < .001
#>   z1                2.993  [ 2.033,  4.407]    < .001
#> 
#> Class 2 ON
#>   Intercept         1.206  [ 0.755,  1.926]     0.432
#>   z1                0.473  [ 0.312,  0.717]    < .001
#> =========================================================
```

### Distal outcomes

A **distal outcome** is a variable that is *caused by* (rather than
constitutive of) the latent classes — for example, a health outcome or a
behavioural measure assessed after class membership has been
established. The 3-step approach first fixes the class solution from the
indicators, then links it to the outcome, preventing the outcome from
distorting the class structure.

**Which bias correction?**

- **BCH** is recommended for continuous distal outcomes. It re-weights
  observations using the inverse of the classification-error matrix,
  without making distributional assumptions about the outcome (Bakk &
  Vermunt, 2016).
- **ML** is recommended for categorical distal outcomes and for
  covariate models. It iterates jointly over the measurement and
  structural components and handles missing outcome data (Vermunt,
  2010).

> **Column convention for `Y`:** column 1 is always the outcome
> variable; any additional columns are treated as covariates.

#### Continuous distal: Outcome ~ Class

Use `structural = "distal_continuous"` to estimate the mean of a
continuous outcome within each class, with no covariate adjustment.
`summary()` reports class-specific means with 95% CIs and robust
sandwich standard errors.

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
#> 
#> Omnibus test (class differences): Wald χ²(2) = 203.09, p  < .001
#> 
#>                  Mean       [95% CI]        SE
#>   Class 1       -0.002  [-0.159,  0.154]     0.080
#>   Class 2        2.875  [ 2.507,  3.244]     0.188
#>   Class 3       -0.782  [-1.430, -0.134]     0.331
#> =========================================================
```

#### Continuous distal with covariate: pooled slope

`structural = "distal_continuous_pooled"` estimates **class-varying
intercepts** (mean differences between classes) but a **single covariate
slope shared across all classes**. This is a parsimonious choice when
you expect the covariate to shift the outcome uniformly, regardless of
class membership.

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
#> Omnibus test (class differences (at covariate zero)): Wald χ²(2) = 321.37, p  < .001
#> 
#>   Latent Class (Intercepts):
#>                  Estimate   [95% CI]        P-Value
#>     Class 1       -0.020  [-0.219,  0.179]     0.845
#>     Class 2        2.895  [ 2.612,  3.178]    < .001
#>     Class 3       -0.778  [-1.122, -0.434]    < .001
#> 
#>   Covariates (Pooled Slopes):
#>                  Estimate   [95% CI]        P-Value
#>     z1            0.018  [-0.108,  0.145]     0.778
#> =========================================================
```

#### Continuous distal with covariate: class-specific slopes

`structural = "distal_continuous_regression"` allows the covariate’s
effect on the outcome to **differ across classes** — i.e., the covariate
is a moderator of the class–outcome relationship. Each class gets its
own intercept and slope.

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
#> Omnibus test (class differences (at covariate zero)): Wald χ²(2) = 385.20, p  < .001
#> 
#> Class 1:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept      -0.040  [-0.240,  0.161]     0.697
#>   z1              0.039  [-0.106,  0.184]     0.599
#> 
#> Class 2:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept       3.545  [ 3.187,  3.903]    < .001
#>   z1              0.628  [ 0.370,  0.887]    < .001
#> 
#> Class 3:
#>                  Estimate   [95% CI]        P-Value
#>   Intercept      -0.984  [-1.307, -0.662]    < .001
#>   z1             -0.865  [-1.161, -0.569]    < .001
#> 
#> ---------------------------------------------------------
#> Wald tests (equality of slopes across classes):
#>                   Wald(χ²(2))  P-Value
#>   z1                 55.66            < .001
#> =========================================================
```

#### Categorical distal: pooled slope

For a binary or ordinal distal outcome, `structural = "distal_pooled"`
estimates class-varying intercepts with a **single covariate slope
shared across all classes**.

``` r
# Simulate a binary categorical distal outcome
set.seed(33)
log_odds <- ifelse(classes == 1, -1.0,
            ifelse(classes == 2,  1.5, 0.0)) + 0.4 * Z[, 1]
D_cat <- matrix(rbinom(n, 1, plogis(log_odds)), ncol = 1)

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
# Reported as odds ratios; class main effects are contrasted against ref_class
summary(fit_pooled)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL DISTAL OUTCOME (POOLED SLOPES)
#> ---------------------------------------------------------
#> 
#> Omnibus test (class differences): Wald χ²(2) = 25.86, p  < .001
#> 
#> Predicted Probabilities (covariates held at zero):
#>                Cat 1    Cat 2   
#>   Class 1        0.695    0.305 
#>   Class 2        0.229    0.771 
#>   Class 3        0.422    0.578 
#> 
#> Pairwise Odds Ratios (Reference: Class 1)
#>                      OR       [95% CI]        P-Value
#> 
#> Outcome Category 2 (vs Category 1) ON
#>   Latent Class:
#>     Class 2          7.683  [ 3.410, 17.310]    < .001
#>     Class 3          3.117  [ 1.403,  6.922]     0.005
#>   Covariates (Pooled Slope):
#>     Z1              1.393  [ 1.109,  1.750]     0.004
#> =========================================================
```

#### Categorical distal: class-specific slopes

`structural = "distal_regression"` relaxes the pooling assumption so
that both the intercept **and** the covariate slope are free to vary by
class. Use this when you have reason to believe the covariate moderates
the class–outcome relationship.

``` r
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
#> CATEGORICAL DISTAL OUTCOME (CLASS-SPECIFIC SLOPES)
#> ---------------------------------------------------------
#> 
#> Predicted Probabilities (covariates held at zero):
#>                Cat 1    Cat 2   
#>   Class 1        0.689    0.311 
#>   Class 2        0.210    0.790 
#>   Class 3        0.424    0.576 
#> 
#> Class-Specific Estimates
#>                      OR       [95% CI]        P-Value
#> 
#> Class 1:
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       0.451  [ 0.238,  0.857]     0.015
#>     Z1              1.341  [ 0.854,  2.105]     0.202
#> 
#> Class 2:
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       3.751  [ 2.579,  5.456]    < .001
#>     Z1              1.570  [ 1.202,  2.052]    < .001
#> 
#> Class 3:
#>   Outcome Category 2 (vs Category 1) ON
#>     Intercept       1.360  [ 0.703,  2.631]     0.361
#>     Z1              1.393  [ 0.781,  2.485]     0.262
#> =========================================================
```

Comparing the pooled and moderated models on AIC / BIC (accessible via
`fit$metrics`) helps determine whether class-specific slopes
meaningfully improve fit:

``` r
cat("Pooled:    BIC =", round(fit_pooled$metrics$bic, 2), "\n")
#> Pooled:    BIC = 703.67
cat("Moderated: BIC =", round(fit_moderated$metrics$bic, 2), "\n")
#> Moderated: BIC = 712.56
```

### Bootstrap Likelihood Ratio Test (BLRT)

The BLRT tests whether a K-class model fits significantly better than a
(K−1)-class model. Because the standard chi-squared approximation is
invalid for mixture models (the null hypothesis places a parameter on
the boundary of the parameter space), the BLRT uses a parametric
bootstrap to build the reference distribution from scratch.

The function returns a list; extract `$p_value` for inference or
`$obs_diff` for the observed test statistic. A histogram of `$null_dist`
with a vertical line at `$obs_diff` provides a visual check.

``` r
result <- calc_blrt(X,
                    k_small     = 2,
                    k_large     = 3,
                    measurement = "binary",
                    n_reps      = 100)
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

# Bootstrap p-value for the 2 vs. 3 class comparison
cat(sprintf("BLRT p-value (2 vs. 3 classes): %.3f\n", result$p_value))
#> BLRT p-value (2 vs. 3 classes): 0.010

# Observed likelihood ratio statistic
cat(sprintf("Observed LR statistic: %.2f\n", result$obs_diff))
#> Observed LR statistic: 54.95
```

## Function reference

| Function | Description |
|----|----|
| `fit_mixture()` | Fit an LCA / LPA model (core function) |
| `compare_mixtures()` | Fit models across a range of K and compare fit indices |
| `print()` | Print a compact model summary |
| `summary()` | Print structural model parameters (OR, CI, p-values; or distal means / regression coefficients) |
| `measurement_summary()` | Print item-response probabilities or profile means per class |
| `classification_diagnostics()` | Average posterior probability (AvePP) matrix |
| `coef()` | Extract odds ratios from a covariate model |
| `confint()` | Confidence intervals for odds ratios (analytical or bootstrap) |
| `analytical_wald_test()` | Analytical Wald chi-squared test for a covariate |
| `bootstrap_covariates()` | Bootstrap standard errors for covariate parameters |
| `wald_omnibus_test()` | Bootstrap Wald omnibus test for a covariate |
| `calc_blrt()` | Bootstrap Likelihood Ratio Test for class enumeration |

## Supported model types

### Measurement (`measurement =`)

| String | Indicator type |
|----|----|
| `"binary"` / `"bernoulli"` | Binary (0/1) |
| `"binary_nan"` / `"bernoulli_nan"` | Binary with missing data |
| `"categorical"` / `"multinoulli"` | Ordinal / polytomous |
| `"categorical_nan"` / `"multinoulli_nan"` | Ordinal with missing data |
| `"continuous"` / `"gaussian_diag"` | Continuous (estimated variance) |
| `"continuous_nan"` / `"gaussian_diag_nan"` | Continuous with missing data |
| `"gaussian"` / `"gaussian_unit"` | Continuous (unit variance) |
| Named list | Mixed model combining any of the above |

### Structural (`structural =`)

| String | Outcome type | Covariate slope | `Y` layout | Recommended correction |
|----|----|----|----|----|
| `"covariate"` | — (predicts class) | — | Covariates only | `"ML"` |
| `"distal_continuous"` | Continuous | None | Outcome only | `"BCH"` |
| `"distal_continuous_pooled"` | Continuous | **Pooled** across classes | Outcome (col 1), covariates (cols 2+) | `"BCH"` |
| `"distal_continuous_regression"` | Continuous | Class-specific | Outcome (col 1), covariates (cols 2+) | `"BCH"` |
| `"distal_pooled"` | Categorical | **Pooled** across classes | Outcome (col 1), covariates (cols 2+) | `"ML"` |
| `"distal_regression"` | Categorical | Class-specific | Outcome (col 1), covariates (cols 2+) | `"ML"` |

> **Pooled vs. moderated:** both continuous and categorical distal
> models come in a pooled (main-effects-only) and a moderated
> (class-specific slopes) version. Start with the pooled model for
> parsimony; switch to the moderated version if you have theoretical
> reasons to expect the covariate effect to differ across classes, and
> verify the improvement using AIC / BIC.

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
    R package. https://github.com/pdvalencia/mixtureEM
