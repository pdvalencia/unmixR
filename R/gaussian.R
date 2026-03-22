# ==============================================================================
# S3 Gaussian Models (Continuous Data)
# ==============================================================================

#' Constructor for Gaussian models
gaussian_model <- function(n_components, type = "gaussian_unit") {
  state <- list(
    n_components = n_components,
    parameters = list()
  )
  class(state) <- c(type, "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# 1. Gaussian Unit (Fixed Variance = 1) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_unit <- function(model_state, X, resp, random_state = NULL) {
  if (!is.null(random_state)) set.seed(random_state)
  idx <- sample.int(nrow(X), model_state$n_components)
  model_state$parameters$means <- X[idx, , drop = FALSE]
  return(model_state)
}

#' @exportS3Method
m_step.gaussian_unit <- function(model_state, X, resp, weights = NULL) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  means <- t(resp) %*% X
  sum_resp <- colSums(resp)
  model_state$parameters$means <- sweep(means, 1, sum_resp, "/")
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_unit <- function(model_state, X) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    ll_matrix <- dnorm(X, mean = mean_c, sd = 1, log = TRUE)
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.gaussian_unit <- function(model_state) {
  return(length(model_state$parameters$means))
}

# ------------------------------------------------------------------------------
# 2. Gaussian Diag (Estimated feature-specific variance) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_diag <- function(model_state, X, resp, random_state = NULL) {
  model_state <- init_params.gaussian_unit(model_state, X, resp, random_state)
  model_state$parameters$covariances <- matrix(1, nrow = model_state$n_components, ncol = ncol(X))
  return(model_state)
}

#' @exportS3Method
m_step.gaussian_diag <- function(model_state, X, resp, weights = NULL) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  # 1. Update Means
  means <- t(resp) %*% X
  sum_resp <- colSums(resp)
  means <- sweep(means, 1, sum_resp, "/")
  model_state$parameters$means <- means

  # 2. Update Variances
  covariances <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))
  for (c in seq_len(model_state$n_components)) {
    diff_sq <- sweep(X, 2, means[c, ], "-")^2
    covariances[c, ] <- colSums(resp[, c] * diff_sq) / sum_resp[c]
  }
  # Regularize to prevent variance collapse
  model_state$parameters$covariances <- covariances + 1e-6
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag <- function(model_state, X) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    sd_c <- matrix(sqrt(model_state$parameters$covariances[c, ]), nrow = n, ncol = ncol(X), byrow = TRUE)

    ll_matrix <- dnorm(X, mean = mean_c, sd = sd_c, log = TRUE)
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.gaussian_diag <- function(model_state) {
  return(length(model_state$parameters$means) + length(model_state$parameters$covariances))
}

# ------------------------------------------------------------------------------
# 3. Gaussian Diag NaN (Missing Data FIML) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method init_params gaussian_diag_nan
init_params.gaussian_diag_nan <- init_params.gaussian_diag
#' @exportS3Method n_parameters gaussian_diag_nan
n_parameters.gaussian_diag_nan <- n_parameters.gaussian_diag

#' @exportS3Method
m_step.gaussian_diag_nan <- function(model_state, X, resp, weights = NULL) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  means <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))
  covariances <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))

  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop=FALSE]
      sum_resp <- colSums(resp_valid)

      # Means
      means[, j] <- t(resp_valid) %*% X[valid, j] / sum_resp

      # Variances
      for (c in seq_len(model_state$n_components)) {
        diff_sq <- (X[valid, j] - means[c, j])^2
        covariances[c, j] <- sum(resp_valid[, c] * diff_sq) / sum_resp[c]
      }
    }
  }
  model_state$parameters$means <- means
  model_state$parameters$covariances <- covariances + 1e-6
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag_nan <- function(model_state, X) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    sd_c <- matrix(sqrt(model_state$parameters$covariances[c, ]), nrow = n, ncol = ncol(X), byrow = TRUE)

    ll_matrix <- dnorm(X, mean = mean_c, sd = sd_c, log = TRUE)
    ll_matrix[is.na(ll_matrix)] <- 0 # FIML Masking
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}
