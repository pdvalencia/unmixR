# ==============================================================================
# S3 Multiparameter Wald Omnibus Test
# ==============================================================================

#' Bootstrap Wald Omnibus Test for a Covariate
#'
#' @description
#' Performs a multiparameter Wald chi-squared test for the effect of a single
#' covariate on latent class membership, using bootstrapped standard errors
#' from \code{\link{bootstrap_covariates}}. Tests the joint null hypothesis
#' that the covariate has no effect across all non-reference classes
#' simultaneously.
#'
#' Compared to \code{\link{analytical_wald_test}}, this test does not rely on
#' Hessian-based standard errors and is therefore more robust in small samples
#' or when classes overlap substantially.
#'
#' @param boot_results A list returned by \code{\link{bootstrap_covariates}}.
#' @param term_name Character string. The name of the covariate to test.
#'   Matched as a substring against the column names of
#'   \code{boot_results$orig_betas} (i.e., partial matches work).
#' @param assume_independence Logical. If \code{TRUE} (default), the Wald
#'   statistic is computed using only the diagonal of the bootstrap covariance
#'   matrix, treating parameters as independent. If \code{FALSE}, the full
#'   covariance matrix is used via the Moore-Penrose pseudoinverse, which
#'   captures correlations between parameters but requires more replications
#'   for stable estimates.
#'
#' @return A single-row data frame with columns `Covariate`, `Wald_Chi2`,
#'   `df`, and `p_value`.
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
#' boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
#' wald_omnibus_test(boot, term_name = "age")
#'
#' # Using the full covariance matrix (requires more bootstrap reps)
#' wald_omnibus_test(boot, term_name = "age", assume_independence = FALSE)
#' }
#'
#' @export
wald_omnibus_test <- function(boot_results, term_name, assume_independence = TRUE) {

  col_names <- colnames(boot_results$orig_betas)
  target_cols <- grep(term_name, col_names)

  if (length(target_cols) == 0) stop("Term not found in model parameters.")

  K <- nrow(boot_results$orig_betas)

  row_sums <- rowSums(abs(boot_results$orig_betas[, target_cols, drop=FALSE]))
  ref_idx <- which(row_sums == 0)[1]
  test_classes <- setdiff(1:K, ref_idx)

  beta_matrix <- boot_results$orig_betas[test_classes, target_cols, drop = FALSE]
  beta_hat <- as.vector(beta_matrix)

  n_reps <- dim(boot_results$boot_betas)[1]
  boot_mat <- matrix(NA, nrow = n_reps, ncol = length(beta_hat))

  col_idx <- 1
  for (v in target_cols) {
    for (c in test_classes) {
      boot_mat[, col_idx] <- boot_results$boot_betas[, c, v]
      col_idx <- col_idx + 1
    }
  }

  Sigma <- cov(boot_mat)

  if (assume_independence) {
    # diag(scalar) in R creates an identity matrix of that size, not a 1x1 matrix.
    # When K=2 there is only one test class so Sigma is 1x1 and diag(Sigma) is a
    # plain scalar. Passing nrow explicitly forces the correct dimensions always.
    inv_Sigma <- diag(1 / diag(Sigma), nrow = length(beta_hat))
  } else {
    inv_Sigma <- pinv(Sigma)
  }

  W_stat <- as.numeric(t(beta_hat) %*% inv_Sigma %*% beta_hat)
  df <- length(beta_hat)
  p_val <- pchisq(W_stat, df = df, lower.tail = FALSE)

  return(data.frame(Covariate = term_name, Wald_Chi2 = round(W_stat, 3), df = df, p_value = round(p_val, 4)))
}
