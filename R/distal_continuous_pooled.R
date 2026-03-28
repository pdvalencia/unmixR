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

  Y_v    <- Y[valid]
  Z_v    <- Z[valid, , drop = FALSE]
  resp_v <- resp[valid, , drop = FALSE]

  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  N_v   <- nrow(Z_v)
  K     <- model_state$n_components
  D_cov <- ncol(Z_v)
  L     <- K + D_cov

  # 1. Expand design matrix for simultaneous intercept/slope estimation
  U <- matrix(0, nrow = N_v * K, ncol = L)
  for (k in 1:K) {
    idx <- ((k - 1) * N_v + 1):(k * N_v)
    U[idx, k] <- 1
    if (D_cov > 0) U[idx, (K + 1):L] <- Z_v
  }

  W_flat <- as.vector(resp_v)
  Y_flat <- rep(Y_v, K)

  # 2. Estimate intercepts and pooled slopes  (bread of the sandwich)
  UWU <- t(U) %*% sweep(U, 1, W_flat, "*")
  UWY <- t(U) %*% (W_flat * Y_flat)

  diag(UWU) <- diag(UWU) + 1e-6        # ridge penalty for stability
  B_inv <- pinv(UWU)
  theta <- B_inv %*% UWY

  # 3. Pooled residual variance
  #
  #    FIX: use SIGNED BCH weights, not abs(W_k).
  #
  #    sigma^2 = sum_k sum_i [ w_ik * (y_i - yhat_ik)^2 ]
  #              / sum_k sum_i w_ik
  #
  #    This matches LatentGOLD's "error variance" under the homoskedastic
  #    (pooled) model and gives a single scalar shared across all classes.
  #    Using abs() inflates sigma^2 and, consequently, all SEs.
  preds  <- U %*% theta
  resids <- as.vector(Y_flat - preds)

  N_total <- sum(W_flat)
  sigma2  <- if (abs(N_total) > 1e-5)
    sum(W_flat * resids^2) / N_total
  else
    1e-5
  sigma2 <- max(sigma2, 1e-5)

  # Store a K x 1 matrix of the pooled variance (same value for every class)
  # to keep the log_likelihood method working unchanged.
  vars <- matrix(sigma2, nrow = K, ncol = 1)

  # 4. Model-based SEs: sqrt( sigma^2 * diag(B^{-1}) )
  #
  #    FIX: replace the naive sandwich   sqrt(diag(B^{-1} M B^{-1}))
  #    with the model-based estimator    sqrt(sigma^2 * diag(B^{-1})).
  #
  #    Rationale: the BCH expanded dataset has K weighted records per
  #    person.  The naive sandwich treats those K records as independent,
  #    inflating the meat by roughly K.  The model-based SE assumes the
  #    normal linear model Y ~ N(X theta, sigma^2 I) with BCH weights,
  #    giving Var(theta) = sigma^2 (U^T W U)^{-1} = sigma^2 B^{-1}.
  #    This reproduces LatentGOLD's reported SEs exactly.
  #
  #    Note: LG reports SEs for CONTRASTS (class k vs. reference class 1)
  #    for the intercepts, while we report SEs for the absolute intercepts.
  #    Both are correct; they differ only in parameterisation.  The
  #    contrast SE is recoverable as sqrt(V[k,k] + V[1,1] - 2*V[k,1])
  #    where V = sigma^2 * B^{-1}.
  ses <- sqrt(pmax(sigma2 * diag(B_inv), 1e-8))

  model_state$parameters$beta_pooled <- matrix(as.vector(theta), nrow = 1)
  model_state$parameters$covariances <- vars
  model_state$parameters$ses         <- matrix(ses, nrow = 1)

  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous_pooled <- function(model_state, X, ...) {
  Y <- as.numeric(X[, 1])
  Z <- impute_covariates(X[, -1, drop = FALSE])

  K     <- model_state$n_components
  D_cov <- ncol(Z)
  L     <- K + D_cov
  N     <- length(Y)
  valid <- !is.na(Y)

  ll    <- matrix(0, nrow = N, ncol = K)
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
                            sd   = sqrt(model_state$parameters$covariances[k, 1]),
                            log  = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous_pooled <- function(model_state, ...) {
  K     <- model_state$n_components
  D_cov <- ncol(model_state$parameters$beta_pooled) - K
  # K intercepts + D_cov slopes + 1 pooled variance
  return(K + D_cov + 1L)
}
