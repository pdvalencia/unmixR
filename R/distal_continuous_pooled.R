# ==============================================================================
# S3 Distal Continuous Pooled Model (Pooled Slopes, Class-Varying Intercepts)
# ==============================================================================

distal_continuous_pooled_model <- function(n_components, ...) {
  state <- list(n_components = n_components, parameters = list())
  class(state) <- c("distal_continuous_pooled", "emission")
  return(state)
}

#' @exportS3Method
init_params.distal_continuous_pooled <- function(model_state, X, resp, ...) {
  Y <- as.numeric(X[, 1])
  Z <- impute_covariates(X[, -1, drop = FALSE])

  K <- model_state$n_components
  D_cov <- ncol(Z)
  L <- K + D_cov

  model_state$parameters$beta_pooled <- matrix(0, nrow = 1, ncol = L)
  # Initialize the intercepts to the global mean
  model_state$parameters$beta_pooled[1, 1:K] <- mean(Y, na.rm = TRUE)

  model_state$parameters$covariances <- matrix(var(Y, na.rm = TRUE), nrow = K, ncol = 1)
  model_state$parameters$ses <- matrix(0, nrow = 1, ncol = L)
  return(model_state)
}

#' @exportS3Method
m_step.distal_continuous_pooled <- function(model_state, X, resp, weights = NULL, ...) {
  Y <- as.numeric(X[, 1])
  Z <- impute_covariates(X[, -1, drop = FALSE])
  valid <- !is.na(Y)

  Y_v <- Y[valid]
  Z_v <- Z[valid, , drop = FALSE]
  resp_v <- resp[valid, , drop = FALSE]

  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  N_v <- nrow(Z_v)
  K <- model_state$n_components
  D_cov <- ncol(Z_v)
  L <- K + D_cov

  # 1. Expand Design Matrix for simultaneous intercept/slope estimation
  U <- matrix(0, nrow = N_v * K, ncol = L)
  for (k in 1:K) {
    idx <- ((k - 1) * N_v + 1):(k * N_v)
    U[idx, k] <- 1
    if (D_cov > 0) U[idx, (K + 1):L] <- Z_v
  }

  W_flat <- as.vector(resp_v)
  Y_flat <- rep(Y_v, K)

  # 2. Estimate Intercepts and Pooled Slopes (The Bread)
  UWU <- t(U) %*% sweep(U, 1, W_flat, "*")
  UWY <- t(U) %*% (W_flat * Y_flat)

  diag(UWU) <- diag(UWU) + 1e-6 # Ridge penalty for stability
  B_inv <- pinv(UWU)
  theta <- B_inv %*% UWY

  # 3. Calculate Class-Specific Variances
  preds <- U %*% theta
  resids <- as.vector(Y_flat - preds)
  vars <- matrix(0, nrow = K, ncol = 1)

  for(k in 1:K) {
    idx <- ((k - 1) * N_v + 1):(k * N_v)
    W_k <- W_flat[idx]
    res_k <- resids[idx]

    Nk_abs <- sum(abs(W_k))
    if (Nk_abs > 1e-5) {
      var_k <- sum(abs(W_k) * res_k^2) / Nk_abs
    } else {
      var_k <- 1e-5
    }
    vars[k, 1] <- max(var_k, 1e-5)
  }

  # 4. Robust Sandwich SEs (The Meat)
  meat_weights <- (W_flat * resids)^2
  M_meat <- t(U) %*% sweep(U, 1, meat_weights, "*")
  cov_matrix <- B_inv %*% M_meat %*% B_inv
  ses <- sqrt(pmax(diag(cov_matrix), 1e-8))

  model_state$parameters$beta_pooled <- matrix(as.vector(theta), nrow = 1)
  model_state$parameters$covariances <- vars
  model_state$parameters$ses <- matrix(ses, nrow = 1)

  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous_pooled <- function(model_state, X, ...) {
  Y <- as.numeric(X[, 1])
  Z <- impute_covariates(X[, -1, drop = FALSE])

  K <- model_state$n_components
  D_cov <- ncol(Z)
  L <- K + D_cov
  N <- length(Y)
  valid <- !is.na(Y)

  ll <- matrix(0, nrow = N, ncol = K)
  theta <- as.vector(model_state$parameters$beta_pooled)

  for (k in 1:K) {
    if (any(valid)) {
      intercept_k <- theta[k]
      if (D_cov > 0) {
        preds <- intercept_k + Z[valid, , drop = FALSE] %*% theta[(K + 1):L]
      } else {
        preds <- rep(intercept_k, sum(valid))
      }

      ll[valid, k] <- dnorm(Y[valid],
                            mean = preds,
                            sd = sqrt(model_state$parameters$covariances[k, 1]),
                            log = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous_pooled <- function(model_state, ...) {
  K <- model_state$n_components
  D_cov <- ncol(model_state$parameters$beta_pooled) - K
  return(K + D_cov + K) # K intercepts + D slopes + K variances
}
