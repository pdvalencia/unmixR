# ==============================================================================
# S3 Distal Continuous Regression (Y ~ Z * Class) - Model-Based SEs
# ==============================================================================

distal_continuous_regression_model <- function(n_components, ...) {
  state <- list(n_components = n_components, parameters = list())
  class(state) <- c("distal_continuous_regression", "emission")
  return(state)
}

#' @exportS3Method
init_params.distal_continuous_regression <- function(model_state, X, resp, ...) {
  Y <- as.numeric(X[, 1])
  Z <- cbind(1, X[, -1, drop = FALSE])
  K <- model_state$n_components
  D <- ncol(Z)
  model_state$parameters$betas       <- matrix(0, nrow = K, ncol = D)
  model_state$parameters$covariances <- matrix(var(Y, na.rm = TRUE), nrow = K, ncol = 1)
  model_state$parameters$ses         <- matrix(0, nrow = K, ncol = D)
  return(model_state)
}

#' @exportS3Method
m_step.distal_continuous_regression <- function(model_state, X, resp,
                                                weights = NULL, ...) {
  Y     <- as.numeric(X[, 1])
  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z     <- cbind(1, Z_raw)
  valid <- !is.na(Y)

  Y_v    <- Y[valid]
  Z_v    <- Z[valid, , drop = FALSE]
  resp_v <- resp[valid, , drop = FALSE]
  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  K <- model_state$n_components
  D <- ncol(Z_v)

  betas  <- matrix(0, nrow = K, ncol = D)
  B_invs <- vector("list", K)      # store for SE computation

  # ── Pass 1: point estimates ────────────────────────────────────────────────
  for (k in seq_len(K)) {
    W_k   <- resp_v[, k]
    ZWZ   <- t(Z_v) %*% sweep(Z_v, 1, W_k, "*")
    ZWY   <- t(Z_v) %*% (W_k * Y_v)
    diag(ZWZ) <- diag(ZWZ) + 1e-6
    B_inv <- pinv(ZWZ)
    betas[k, ]  <- as.vector(B_inv %*% ZWY)
    B_invs[[k]] <- B_inv
  }

  # ── Pass 2: pooled residual variance (SIGNED BCH weights) ─────────────────
  #
  #    FIX: use signed weights, not abs(W_k), and pool across ALL classes.
  #
  #    sigma^2 = sum_k sum_i [ w_ik * (y_i - yhat_ik)^2 ]
  #              / sum_k sum_i w_ik
  #
  #    This matches LatentGOLD's pooled "error variance" (homoskedastic model).
  #    Using abs() inflates sigma^2 and all downstream SEs.
  total_ss <- 0
  total_n  <- 0
  for (k in seq_len(K)) {
    W_k      <- resp_v[, k]
    resids_k <- Y_v - Z_v %*% betas[k, ]
    total_ss <- total_ss + sum(W_k * resids_k^2)
    total_n  <- total_n  + sum(W_k)
  }
  sigma2 <- max(total_ss / total_n, 1e-5)

  # Store as K x 1 matrix (same value repeated) to keep log_likelihood
  # working without changes.
  vars <- matrix(sigma2, nrow = K, ncol = 1)

  # ── Pass 3: model-based SEs ────────────────────────────────────────────────
  #
  #    FIX: replace naive sandwich sqrt(diag(B^{-1} M B^{-1})) with
  #    model-based sqrt(sigma^2 * diag(B_inv_k)).
  #
  #    Each class has its own B_inv_k (separate WLS), so the SE for class k
  #    is sqrt(sigma^2 * diag(B_inv_k)).  The pooled sigma^2 is shared.
  #    This exactly reproduces LatentGOLD's reported SEs.
  ses <- matrix(0, nrow = K, ncol = D)
  for (k in seq_len(K))
    ses[k, ] <- sqrt(pmax(sigma2 * diag(B_invs[[k]]), 1e-8))

  model_state$parameters$betas       <- betas
  model_state$parameters$covariances <- vars
  model_state$parameters$ses         <- ses

  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous_regression <- function(model_state, X, ...) {
  Y     <- as.numeric(X[, 1])
  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z     <- cbind(1, Z_raw)
  K     <- model_state$n_components
  N     <- length(Y)
  ll    <- matrix(0, nrow = N, ncol = K)
  valid <- !is.na(Y)
  for (k in seq_len(K)) {
    if (any(valid)) {
      preds <- Z[valid, , drop = FALSE] %*% model_state$parameters$betas[k, ]
      ll[valid, k] <- dnorm(Y[valid],
                            mean = preds,
                            sd   = sqrt(model_state$parameters$covariances[k, 1]),
                            log  = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous_regression <- function(model_state, ...) {
  K <- model_state$n_components
  D <- ncol(model_state$parameters$betas)
  # K*D regression coefficients + 1 pooled variance
  return(K * D + 1L)
}
