# ==============================================================================
# Analytical Wald Omnibus Test
# ==============================================================================

#' Analytical Wald Test for a Single Covariate
#'
#' @description
#' Performs an omnibus Wald chi-squared test for the effect of a single
#' covariate on latent class membership, using standard errors derived
#' analytically from the observed information matrix (Hessian). The null
#' hypothesis is that the covariate has no effect on class membership across
#' all non-reference classes simultaneously.
#'
#' For small samples or poorly conditioned Hessians, consider
#' \code{\link{wald_omnibus_test}} instead, which uses bootstrapped standard
#' errors from \code{\link{bootstrap_covariates}}.
#'
#' @param model A fitted \code{mixture_model} object with a covariate
#'   structural model (fitted with \code{structural = "covariate"}).
#' @param term_name Character string. The name of the covariate to test.
#'   Matched as a substring against the column names of the model's beta
#'   matrix (i.e., \code{grep} is used internally), so partial matches work.
#' @param ref_class Integer. The reference latent class. Default is \code{1}.
#'
#' @return A single-row data frame with columns `Covariate`, `Wald_Chi2`,
#'   `df`, `p_value`, and `Method` (always `"Analytical"`).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 3, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' analytical_wald_test(fit, term_name = "age")
#' }
#'
#' @export
analytical_wald_test <- function(model, term_name, ref_class = 1) {
  if (is.null(model$sm) || !inherits(model$sm, "covariate")) stop("No covariate model.")

  K <- model$n_components

  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  # K=1 means no contrasts can be formed (there is only one class).
  if (K == 1)
    stop("analytical_wald_test requires at least 2 classes. Got n_components = 1.")

  H <- model$sm$parameters$hessian
  if (is.null(H) || all(H == 0)) stop("Hessian matrix is missing. Refit the model.")

  Sigma_full <- pinv(-H)
  K <- model$n_components
  D <- ncol(model$sm$parameters$beta)
  col_names <- colnames(model$sm$parameters$beta)
  target_col_idx <- grep(term_name, col_names)

  if (length(target_col_idx) == 0) stop("Variable name not found.")

  total_params <- K * D
  test_classes <- setdiff(1:K, ref_class)
  num_tests <- length(test_classes) * length(target_col_idx)

  C <- matrix(0, nrow = num_tests, ncol = total_params)
  row <- 1
  for (c in test_classes) {
    for (v in target_col_idx) {
      idx_current <- (c - 1) * D + v
      idx_ref <- (ref_class - 1) * D + v
      C[row, idx_current] <- 1
      C[row, idx_ref] <- -1
      row <- row + 1
    }
  }

  beta_vec <- as.vector(t(model$sm$parameters$beta))
  diff_beta <- C %*% beta_vec
  V_cov <- C %*% Sigma_full %*% t(C)

  # Safe matrix inversion
  W_stat <- as.numeric(t(diff_beta) %*% pinv(V_cov) %*% diff_beta)
  df <- num_tests
  p_val <- pchisq(W_stat, df = df, lower.tail = FALSE)

  result <- data.frame(
    Covariate = term_name,
    Wald_Chi2 = round(W_stat, 3),
    df        = df,
    p_value   = p_val,
    Method    = "Analytical"
  )
  class(result) <- c("mixture_wald", "data.frame")
  return(result)
}

#' @export
print.mixture_wald <- function(x, ...) {
  cat("=========================================================\n")
  cat("                 WALD TEST (COVARIATE)                   \n")
  cat("=========================================================\n")
  cat(sprintf("  Covariate : %s\n", x$Covariate))
  cat(sprintf("  Method    : %s\n", x$Method))
  cat("---------------------------------------------------------\n")
  cat(sprintf("  Wald \u03c7\u00b2(%d) = %.3f,  p%s\n",
              x$df, x$Wald_Chi2,
              if (x$p_value < 0.001) " < .001"
              else sprintf(" = %.3f", x$p_value)))
  cat("=========================================================\n")
  invisible(x)
}
