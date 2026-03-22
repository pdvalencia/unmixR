#' Confidence Intervals for Odds Ratios in a Mixture Model
#'
#' @description
#' Computes confidence intervals for the odds ratios of covariate effects on
#' latent class membership. By default, analytical standard errors derived
#' from the observed information matrix (Hessian) are used. Bootstrapped
#' standard errors from \code{\link{bootstrap_covariates}} can be supplied
#' instead, which is recommended in small samples or when the Hessian is
#' poorly conditioned.
#'
#' @param object A fitted \code{mixture_model} object with a covariate
#'   structural model (fitted with \code{structural = "covariate"}).
#' @param boot_results Optional. A list returned by
#'   \code{\link{bootstrap_covariates}}. When provided, bootstrapped standard
#'   errors are used. When \code{NULL} (default), analytical SEs derived from
#'   the Hessian are used.
#' @param ref_class Integer. The reference latent class. Odds ratios for all
#'   other classes are expressed relative to this class. Default is \code{1}.
#' @param level Numeric. The confidence level. Default is \code{0.95}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return A named list with one data frame per predictor variable (including
#'   the intercept). Each data frame has columns `OR` (odds ratio), `Lower`,
#'   and `Upper` (confidence bounds), with one row per latent class. Values
#'   are rounded to three decimal places.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#'
#' # Analytical CIs (from Hessian)
#' confint(fit)
#'
#' \dontrun{
#' # Bootstrapped CIs
#' boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
#' confint(fit, boot_results = boot)
#' }
#'
#' @export
confint.mixture_model <- function(object, boot_results = NULL, ref_class = 1, level = 0.95, ...) {
  if (is.null(object$sm) || !inherits(object$sm, "covariate")) stop("No covariate model.")

  K <- object$n_components

  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  D <- ncol(object$sm$parameters$beta)
  z_crit <- qnorm(1 - (1 - level) / 2)

  # 1. Extract Beta (centered on ref_class)
  betas <- object$sm$parameters$beta
  betas_ref <- sweep(betas, 2, betas[ref_class, ], "-")

  # 2. Determine Standard Errors (Analytical vs Bootstrapped)
  if (!is.null(boot_results)) {
    message("Calculating Bootstrapped Confidence Intervals...")
    se <- boot_results$standard_errors # Already centered by our bootstrap function
  } else {
    message("Calculating Analytical Confidence Intervals (Hessian)...")
    H <- object$sm$parameters$hessian
    if (is.null(H)) stop("Hessian missing. Refit model.")

    # Extract SEs from the diagonal of the inverted Hessian
    # Note: Analytical SEs are for the absolute parameters;
    # for differences (vs ref), we use the contrast matrix logic
    Sigma <- pinv(-H)
    se <- matrix(0, nrow = K, ncol = D)
    for (c in 1:K) {
      for (v in 1:D) {
        # Variance of a difference: Var(B_c - B_ref) = Var(B_c) + Var(B_ref) - 2*Cov(B_c, B_ref)
        idx_c <- (c - 1) * D + v
        idx_ref <- (ref_class - 1) * D + v
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] - 2 * Sigma[idx_c, idx_ref]
        se[c, v] <- sqrt(pmax(0, var_diff))
      }
    }
  }

  # 3. Calculate Bounds
  lower_beta <- betas_ref - z_crit * se
  upper_beta <- betas_ref + z_crit * se

  # 4. Transform to Odds Ratios
  OR <- exp(betas_ref)
  LB <- exp(lower_beta)
  UB <- exp(upper_beta)

  # 5. Format Output
  res_list <- list()
  cov_names <- colnames(betas)
  if (is.null(cov_names)) cov_names <- paste0("V", seq_len(D))

  for (v in 1:D) {
    df <- data.frame(
      OR = OR[, v],
      Lower = LB[, v],
      Upper = UB[, v]
    )
    rownames(df) <- paste("Class", 1:K)
    rownames(df)[ref_class] <- paste("Class", ref_class, "(Ref)")
    res_list[[cov_names[v]]] <- round(df, 3)
  }

  return(res_list)
}
