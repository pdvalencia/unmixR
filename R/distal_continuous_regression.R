# ==============================================================================
# S3 Distal Continuous Regression (Y ~ Z * Class) - ROBUST SEs
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

  model_state$parameters$betas <- matrix(0, nrow = K, ncol = D)
  model_state$parameters$covariances <- matrix(var(Y, na.rm=TRUE), nrow = K, ncol = 1)
  model_state$parameters$ses <- matrix(0, nrow = K, ncol = D)
  return(model_state)
}

#' @exportS3Method
m_step.distal_continuous_regression <- function(model_state, X, resp, weights = NULL) {
  Y <- as.numeric(X[, 1])
  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z <- cbind(1, Z_raw)
  valid <- !is.na(Y)

  Y_v <- Y[valid]
  Z_v <- Z[valid, , drop=FALSE]
  resp_v <- resp[valid, , drop=FALSE]

  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  K <- model_state$n_components
  D <- ncol(Z_v)

  betas <- matrix(0, nrow = K, ncol = D)
  vars <- matrix(0, nrow = K, ncol = 1)
  ses <- matrix(0, nrow = K, ncol = D)

  for (k in 1:K) {
    W_k <- resp_v[, k]

    # --- 1. THE BREAD (Point Estimates) ---
    ZWZ <- t(Z_v) %*% sweep(Z_v, 1, W_k, "*")
    ZWY <- t(Z_v) %*% (W_k * Y_v)

    diag(ZWZ) <- diag(ZWZ) + 1e-6 # Ridge penalty for stability
    B_inv <- pinv(ZWZ)

    beta_k <- B_inv %*% ZWY
    betas[k, ] <- as.vector(beta_k)

    # --- 2. THE RESIDUALS ---
    preds <- Z_v %*% beta_k
    resids <- as.vector(Y_v - preds)

    # Calculate scalar variance for the log-likelihood using ABSOLUTE weights
    # This keeps the shape of the normal distribution valid and prevents NaN crashes
    Nk_abs <- sum(abs(W_k))
    if (Nk_abs > 1e-5) {
      var_k <- sum(abs(W_k) * resids^2) / Nk_abs
    } else {
      var_k <- 1e-5
    }
    vars[k, 1] <- max(var_k, 1e-5)

    # --- 3. THE MEAT (Robust Sandwich SEs) ---
    # Meat = Z' * diag(W^2 * e^2) * Z
    meat_weights <- (W_k * resids)^2
    M_meat <- t(Z_v) %*% sweep(Z_v, 1, meat_weights, "*")

    # Sandwich: Cov(Beta) = B_inv * M_meat * B_inv
    cov_matrix <- B_inv %*% M_meat %*% B_inv
    ses[k, ] <- sqrt(pmax(diag(cov_matrix), 1e-8)) # Protect against floating point zeros
  }

  model_state$parameters$betas <- betas
  model_state$parameters$covariances <- vars
  model_state$parameters$ses <- ses
  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous_regression <- function(model_state, X) {
  Y <- as.numeric(X[, 1])
  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z <- cbind(1, Z_raw)

  K <- model_state$n_components
  N <- length(Y)
  ll <- matrix(0, nrow = N, ncol = K)
  valid <- !is.na(Y)

  for (k in 1:K) {
    if (any(valid)) {
      beta_k <- model_state$parameters$betas[k, ]
      preds <- Z[valid, , drop=FALSE] %*% beta_k

      ll[valid, k] <- dnorm(Y[valid],
                            mean = preds,
                            sd = sqrt(model_state$parameters$covariances[k, 1]),
                            log = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous_regression <- function(model_state) {
  K <- model_state$n_components
  D <- ncol(model_state$parameters$betas)
  return((K * D) + K)
}
