# ==============================================================================
# Bootstrap Likelihood Ratio Test (BLRT) - DYNAMIC MEASUREMENT
# ==============================================================================

# ------------------------------------------------------------------------------
# HELPER: Dynamically generate synthetic data from a fitted measurement model
# (Uses "Duck Typing" to perfectly handle any S3 class name)
# ------------------------------------------------------------------------------
generate_synthetic_data <- function(mm, classes, N) {
  if (inherits(mm, "nested")) {
    res <- list()
    for (name in names(mm$models)) {
      res[[name]] <- generate_synthetic_data(mm$models[[name]], classes, N)
    }
    return(do.call(cbind, res))
  }

  is_cont <- !is.null(mm$parameters$means)
  is_poly <- !is.null(mm$max_val)
  is_bin  <- !is.null(mm$parameters$pis) && is.null(mm$max_val)

  if (is_cont) {
    D <- ncol(mm$parameters$means)
    X_gen <- matrix(0, nrow = N, ncol = D)
    for (i in 1:N) {
      vars <- if (!is.null(mm$parameters$covariances)) mm$parameters$covariances[classes[i], ] else rep(1, D)
      X_gen[i, ] <- rnorm(D, mean = mm$parameters$means[classes[i], ], sd = sqrt(vars))
    }
    return(X_gen)

  } else if (is_poly) {
    D <- ncol(mm$parameters$pis)
    M <- mm$max_val
    n_items <- D / M
    X_gen <- matrix(0L, nrow = N, ncol = n_items)
    for (i in 1:N) {
      probs <- matrix(mm$parameters$pis[classes[i], ], nrow = n_items, ncol = M, byrow = TRUE)
      for (j in 1:n_items) {
        X_gen[i, j] <- sample(1:M, 1, prob = probs[j, ])
      }
    }
    return(X_gen)

  } else if (is_bin) {
    D <- ncol(mm$parameters$pis)
    X_gen <- matrix(0, nrow = N, ncol = D)
    for (i in 1:N) {
      X_gen[i, ] <- rbinom(D, 1, prob = mm$parameters$pis[classes[i], ])
    }
    return(X_gen)

  } else {
    stop(sprintf("Unsupported measurement model parameters. S3 Class: %s", paste(class(mm), collapse=", ")))
  }
}

# ------------------------------------------------------------------------------
# MAIN BLRT FUNCTION
# ------------------------------------------------------------------------------

#' Bootstrap Likelihood Ratio Test (BLRT) for Class Enumeration
#'
#' @description
#' Tests whether a K-class model fits significantly better than a (K-1)-class
#' model by using a parametric bootstrap to approximate the null distribution
#' of the likelihood ratio statistic. This avoids the known violation of
#' standard chi-squared regularity conditions in mixture models, where the
#' null hypothesis places a parameter on the boundary of the parameter space.
#'
#' The function fits both the null (\code{k_small}-class) and alternative
#' (\code{k_large}-class) models on the observed data, computes the observed
#' likelihood ratio statistic, then generates \code{n_reps} synthetic datasets
#' under the null model to build the reference distribution.
#'
#' @param X A numeric matrix of indicator variables (the measurement data).
#'   Rows are observations; columns are items.
#' @param k_small Positive integer. The number of classes in the null model.
#' @param k_large Positive integer. The number of classes in the alternative
#'   model. Must be strictly greater than \code{k_small}.
#' @param measurement Character string specifying the measurement model type.
#'   See \code{\link{fit_mixture}} for accepted values. Default is
#'   \code{"binary"}.
#' @param n_reps Positive integer. Number of bootstrap replications. More
#'   replications give a more precise p-value but take longer. Default is
#'   \code{100}.
#' @param n_init_base Positive integer. Number of random restarts when fitting
#'   the observed-data models. Default is \code{20}.
#' @param n_init_boot Positive integer. Number of random restarts for each
#'   bootstrap replicate pair of models. Default is \code{10}.
#' @param ... Additional arguments passed to \code{\link{fit_mixture}}.
#'
#' @return A named list with three elements:
#'   * `p_value` Numeric. The bootstrap p-value, computed as
#'     `(#{null_dist >= obs_diff} + 1) / (n_reps + 1)`.
#'   * `obs_diff` Numeric. The observed likelihood ratio statistic:
#'     `2 * (LL_large - LL_small)`.
#'   * `null_dist` Numeric vector of length `n_reps` containing the bootstrap
#'     null distribution of the test statistic.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' result <- calc_blrt(X, k_small = 2, k_large = 3,
#'                     measurement = "binary", n_reps = 50)
#' result$p_value
#' hist(result$null_dist, main = "BLRT Null Distribution")
#' abline(v = result$obs_diff, col = "red")
#' }
#'
#' @export
calc_blrt <- function(X, k_small, k_large, measurement = "binary",
                      n_reps = 100, n_init_base = 20, n_init_boot = 10, ...) {

  # k_large must be strictly greater than k_small.
  # When k_small == k_large the two independently fitted models can land in
  # different local optima, producing a non-zero obs_diff that is meaningless
  # as a likelihood ratio test statistic.
  if (k_small >= k_large)
    stop(sprintf(
      "k_large (%d) must be strictly greater than k_small (%d) for the BLRT.",
      k_large, k_small
    ))

  message(sprintf("Running BLRT: %d vs %d classes (This may take a while)...", k_small, k_large))

  null_model <- fit_mixture(X, n_components = k_small, measurement = measurement, n_init = n_init_base, ...)
  alt_model  <- fit_mixture(X, n_components = k_large, measurement = measurement, n_init = n_init_base, ...)

  obs_diff <- 2 * (alt_model$metrics$ll - null_model$metrics$ll)
  null_dist <- numeric(n_reps)
  N <- nrow(X)

  for (i in 1:n_reps) {
    classes <- sample(1:k_small, size = N, replace = TRUE, prob = null_model$weights)
    X_gen <- generate_synthetic_data(null_model$mm, classes, N)

    # refine = FALSE: bootstrap replicates only need the likelihood ratio, not
    # polished final estimates. Skipping L-BFGS makes each replicate ~100x faster
    # with no effect on the validity of the p-value.
    m_null_gen <- fit_mixture(X_gen, n_components = k_small, measurement = measurement,
                              n_init = n_init_boot, refine = FALSE, ...)
    m_alt_gen  <- fit_mixture(X_gen, n_components = k_large, measurement = measurement,
                              n_init = n_init_boot, refine = FALSE, ...)

    null_dist[i] <- 2 * (m_alt_gen$metrics$ll - m_null_gen$metrics$ll)

    if (i %% 10 == 0) cat(sprintf("  Bootstrap draw %d / %d completed.\n", i, n_reps))
  }

  p_val <- (sum(null_dist >= obs_diff) + 1) / (n_reps + 1)

  return(list(
    p_value = p_val,
    obs_diff = obs_diff,
    null_dist = null_dist
  ))
}
