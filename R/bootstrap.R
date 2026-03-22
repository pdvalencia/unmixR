# ==============================================================================
# S3 Bootstrapping with Optimal Label-Switching Alignment
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal helper: all permutations of 1:n as a list of integer vectors.
# Used by align_classes() for globally optimal assignment.
# ------------------------------------------------------------------------------
.permutations <- function(n) {
  if (n == 1L) return(list(1L))
  sub    <- .permutations(n - 1L)
  result <- vector("list", factorial(n))
  k      <- 1L
  for (ins in 1:n) {
    for (p in sub) {
      result[[k]] <- c(ins, p + (p >= ins))
      k <- k + 1L
    }
  }
  result
}

# ------------------------------------------------------------------------------
# Internal helper: extract a K x P numeric matrix that characterises each
# class and can be used for label-switching alignment.
#
# For Bernoulli/categorical models this is the item-probability matrix ($pis).
# For Gaussian models this is the means matrix ($means).
# For nested models the sub-model matrices are column-bound together so that
# all class-specific information contributes to the alignment.
# ------------------------------------------------------------------------------
get_mm_alignment_matrix <- function(mm) {
  if (inherits(mm, "nested")) {
    # Combine every sub-model's alignment matrix column-wise.
    parts <- lapply(mm$models, get_mm_alignment_matrix)
    return(do.call(cbind, parts))
  }
  if (!is.null(mm$parameters[["pis"]]))   return(mm$parameters$pis)
  if (!is.null(mm$parameters[["means"]])) return(mm$parameters$means)
  stop(paste(
    "Cannot determine alignment matrix: measurement model of class",
    paste(class(mm), collapse = "/"),
    "has neither $pis nor $means parameters."
  ))
}

# Internal helper: align bootstrap classes to original model classes via
# optimal assignment. For K <= 8 the assignment is solved exactly by
# enumerating all K! permutations of the cost matrix. For K > 8 a greedy
# nearest-neighbour fallback is used.
#
# orig_mat  K x P matrix characterising each class in the original model
# boot_mat  K x P matrix characterising each class in a bootstrap model
# Returns an integer vector of length K: mapping[i] = which boot class
# matches orig class i.
align_classes <- function(orig_mat, boot_mat) {
  K <- nrow(orig_mat)

  # Squared-distance cost matrix: cost[i, j] = dist(orig_i, boot_j)
  cost <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K)
    for (j in 1:K)
      cost[i, j] <- sum((orig_mat[i, ] - boot_mat[j, ])^2)

  if (K <= 8) {
    # Globally optimal: enumerate all permutations
    best_cost <- Inf
    best_perm <- 1:K
    for (perm in .permutations(K)) {
      perm    <- unlist(perm)
      c_val   <- sum(cost[cbind(1:K, perm)])
      if (c_val < best_cost) {
        best_cost <- c_val
        best_perm <- perm
      }
    }
    return(best_perm)
  }

  # Greedy fallback for K > 8
  mapping   <- integer(K)
  available <- 1:K
  for (i in 1:K) {
    j          <- which.min(cost[i, available])
    mapping[i] <- available[j]
    available  <- available[-j]
  }
  return(mapping)
}

#' Bootstrap Standard Errors for Covariate Model Parameters
#'
#' @description
#' Estimates standard errors and p-values for the covariate structural model
#' parameters using nonparametric bootstrap resampling. In each replicate,
#' a new model is fitted on a bootstrap sample and class labels are aligned
#' to those of the original model via globally optimal assignment (solving the
#' linear sum assignment problem by enumeration for K <= 8), correcting for
#' label switching. Results can be passed directly to
#' \code{\link{confint.mixture_model}} and \code{\link{wald_omnibus_test}}.
#'
#' @param model_state A fitted \code{mixture_model} object with a covariate
#'   structural model (fitted with \code{structural = "covariate"}).
#' @param X Numeric matrix. The original measurement data used to fit
#'   \code{model_state}.
#' @param Y Numeric matrix. The original covariate data used to fit
#'   \code{model_state}.
#' @param n_reps Positive integer. Number of bootstrap replications.
#'   Default is \code{100}.
#' @param random_state Integer seed for reproducibility. Default is \code{123}.
#' @param ref_class Integer. Reference class for centering bootstrap betas.
#'   Should match the \code{ref_class} used in subsequent calls to
#'   \code{\link{confint.mixture_model}} and \code{\link{wald_omnibus_test}}.
#'   Default is \code{1}.
#'
#' @return A named list with four elements:
#'   * `standard_errors` K x D numeric matrix of bootstrap standard errors
#'     (one row per class, one column per predictor).
#'   * `p_values` K x D numeric matrix of two-sided p-values.
#'   * `boot_betas` 3-D array of shape `(n_reps, K, D)` containing aligned,
#'     reference-centered bootstrap beta estimates. Required by
#'     [wald_omnibus_test()].
#'   * `orig_betas` K x D matrix of the original model's beta coefficients,
#'     centered on `ref_class`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
#' boot$p_values
#' confint(fit, boot_results = boot)
#' }
#'
#' @export
bootstrap_covariates <- function(model_state, X, Y, n_reps = 100,
                                 random_state = 123, ref_class = 1) {
  set.seed(random_state)

  n_samples <- nrow(X)
  K         <- model_state$n_components
  D_cov     <- ncol(model_state$sm$parameters$beta)

  # Extract alignment matrix from the original model's measurement model.
  # get_mm_alignment_matrix() handles Bernoulli ($pis), Gaussian ($means), and
  orig_align_mat <- get_mm_alignment_matrix(model_state$mm)

  # Centre original betas on the reference class
  orig_betas <- sweep(model_state$sm$parameters$beta, 2,
                      model_state$sm$parameters$beta[ref_class, ], "-")
  y_names <- colnames(Y)
  if (is.null(y_names)) y_names <- paste0("V", seq_len(ncol(Y)))
  colnames(orig_betas) <- c("Intercept", y_names)

  boot_betas <- array(NA, dim = c(n_reps, K, D_cov))

  cat(sprintf("Running bootstrap with alignment (Ref Class: %d, %d reps)...\n",
              ref_class, n_reps))
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)

  # Retrieve the original measurement descriptor stored by fit_mixture().
  # Using class(model_state$mm)[1] was wrong for nested models: it returns
  # "nested", which build_emission() does not recognise as a string descriptor.
  # The stored descriptor (a list or a string) faithfully reproduces the
  # original measurement specification in each bootstrap replicate (Bug 2 fix).
  measurement_desc <- model_state$measurement_descriptor
  if (is.null(measurement_desc))
    stop(paste(
      "model_state$measurement_descriptor is NULL.",
      "Refit the original model with the current version of fit_mixture(),",
      "which stores the descriptor automatically."
    ))

  for (i in 1:n_reps) {
    idx    <- sample(1:n_samples, replace = TRUE)
    X_boot <- X[idx, , drop = FALSE]
    Y_boot <- Y[idx, , drop = FALSE]

    b_model <- fit_mixture(
      X = X_boot, Y = Y_boot, n_components = K,
      measurement  = measurement_desc,   # use stored descriptor (Bug 2 fix)
      structural   = "covariate",
      n_steps      = 3,
      correction   = "ML",
      n_init       = 1,
      order_by_size = FALSE   # preserve internal order for alignment
    )

    # Align bootstrap classes to original model using get_mm_alignment_matrix()
    # so Gaussian and nested models are handled correctly (Bug 3 fix).
    boot_align_mat <- get_mm_alignment_matrix(b_model$mm)
    mapping        <- align_classes(orig_align_mat, boot_align_mat)
    aligned_betas  <- b_model$sm$parameters$beta[mapping, ]

    # Centre on the reference class
    aligned_betas_ref <- sweep(aligned_betas, 2,
                               aligned_betas[ref_class, ], "-")
    boot_betas[i, , ] <- aligned_betas_ref
    setTxtProgressBar(pb, i)
  }
  close(pb)

  se_betas <- apply(boot_betas, c(2, 3), sd)
  z_scores <- orig_betas / se_betas
  p_values <- 2 * (1 - pnorm(abs(z_scores)))

  rownames(p_values) <- paste("Class", 1:K)
  rownames(p_values)[ref_class] <- paste("Class", ref_class, "(Ref)")

  return(list(
    standard_errors = se_betas,
    p_values        = p_values,
    boot_betas      = boot_betas,
    orig_betas      = orig_betas
  ))
}
