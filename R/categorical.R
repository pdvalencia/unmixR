# ==============================================================================
# S3 Categorical Models (Bernoulli and Multinoulli)
# ==============================================================================

#' Constructor for Categorical models
categorical_model <- function(n_components, type = "bernoulli", max_val = NULL, ...) {
  state <- list(
    n_components = n_components,
    parameters = list(),
    max_val = max_val
  )
  class(state) <- c(type, "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# 1. Bernoulli (Binary) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.bernoulli <- function(model_state, X, resp, random_state = NULL) {
  if (!is.null(random_state)) set.seed(random_state)
  model_state$parameters$pis <- matrix(
    runif(model_state$n_components * ncol(X), 0.25, 0.75),
    nrow = model_state$n_components
  )
  return(model_state)
}

#' @exportS3Method
m_step.bernoulli <- function(model_state, X, resp, weights = NULL, alpha = 1.0) {
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
    marginal_prob <- colSums(sweep(X, 1, weights, "*"), na.rm = TRUE) / sum(weights)
  } else {
    marginal_prob <- colMeans(X, na.rm = TRUE)
  }

  K <- model_state$n_components
  prior_obs <- alpha / K

  pis <- t(resp) %*% X

  # Add Dirichlet prior to numerator (proportional to marginal probability)
  pis <- sweep(pis, 2, prior_obs * marginal_prob, "+")

  # Add Dirichlet prior to denominator
  sum_resp <- colSums(resp) + prior_obs
  pis <- sweep(pis, 1, sum_resp, "/")

  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.bernoulli <- function(model_state, X) {
  log_eps <- matrix(0, nrow = nrow(X), ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    pi_c <- model_state$parameters$pis[c, ]
    # Small epsilon just in case, though the prior makes 0s impossible
    pi_c <- pmax(pmin(pi_c, 1 - 1e-15), 1e-15)
    term1 <- sweep(X, 2, log(pi_c), "*")
    term2 <- sweep(1 - X, 2, log(1 - pi_c), "*")
    log_eps[, c] <- rowSums(term1 + term2)
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.bernoulli <- function(model_state) {
  return(length(model_state$parameters$pis))
}

# ------------------------------------------------------------------------------
# 2. Bernoulli NaN (Binary with Missing Data) S3 Methods
# ------------------------------------------------------------------------------
#' @exportS3Method init_params bernoulli_nan
init_params.bernoulli_nan <- init_params.bernoulli
#' @exportS3Method n_parameters bernoulli_nan
n_parameters.bernoulli_nan <- n_parameters.bernoulli

#' @exportS3Method
m_step.bernoulli_nan <- function(model_state, X, resp, weights = NULL, alpha = 1.0) {
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
  }

  K <- model_state$n_components
  prior_obs <- alpha / K
  pis <- matrix(0, nrow = K, ncol = ncol(X))

  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop = FALSE]

      if (!is.null(weights)) {
        marg_prob <- sum(X[valid, j] * weights[valid]) / sum(weights[valid])
      } else {
        marg_prob <- mean(X[valid, j])
      }

      num <- t(resp_valid) %*% X[valid, j] + (prior_obs * marg_prob)
      den <- colSums(resp_valid) + prior_obs
      pis[, j] <- num / den
    }
  }
  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.bernoulli_nan <- function(model_state, X) {
  log_eps <- matrix(0, nrow = nrow(X), ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    pi_c <- model_state$parameters$pis[c, ]
    pi_c <- pmax(pmin(pi_c, 1 - 1e-15), 1e-15)
    term1 <- sweep(X, 2, log(pi_c), "*")
    term2 <- sweep(1 - X, 2, log(1 - pi_c), "*")
    sum_mat <- term1 + term2
    sum_mat[is.na(sum_mat)] <- 0 # FIML: Ignore missing values in the sum
    log_eps[, c] <- rowSums(sum_mat)
  }
  return(log_eps)
}

# ------------------------------------------------------------------------------
# 3. Multinoulli (Categorical) S3 Methods
# ------------------------------------------------------------------------------

one_hot <- function(X, max_val) {
  n <- nrow(X)
  D <- ncol(X)
  out <- matrix(0, nrow = n, ncol = D * max_val)

  row_indices <- row(X)
  col_base <- (col(X) - 1) * max_val
  valid <- !is.na(X)

  # Guard against float values: R silently truncates non-integer subscripts,
  # so 1.9 is treated as category 1 and 2.9 as category 2 with no error.
  # Catch this before it produces silent wrong results.
  if (any(valid)) {
    valid_vals <- X[valid]
    if (!all(valid_vals == as.integer(valid_vals)))
      stop(paste(
        "one_hot: X contains non-integer values in a categorical column.",
        "Categorical items must be integer-valued (e.g. 1L, 2L, 3L).",
        "Did you forget to round or convert to integer?"
      ))
  }

  if (any(valid)) {
    target_rows <- row_indices[valid]
    target_cols <- col_base[valid] + X[valid]
    out[cbind(target_rows, target_cols)] <- 1
  }
  return(out)
}

#' @exportS3Method
init_params.multinoulli <- function(model_state, X, resp, random_state = NULL) {
  if (!is.null(random_state)) set.seed(random_state)

  # Only infer max_val from the data if the user didn't explicitly provide it
  if (is.null(model_state$max_val)) {
    model_state$max_val <- max(X, na.rm = TRUE)
  }

  n_features <- ncol(X) * model_state$max_val
  pis <- matrix(runif(model_state$n_components * n_features), nrow = model_state$n_components)

  for (j in seq_len(ncol(X))) {
    cols <- ((j - 1) * model_state$max_val + 1):(j * model_state$max_val)
    pis[, cols] <- sweep(pis[, cols, drop=FALSE], 1, rowSums(pis[, cols, drop=FALSE]), "/")
  }
  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
m_step.multinoulli <- function(model_state, X, resp, weights = NULL, alpha = 1.0) {
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
  }

  X_oh <- one_hot(X, model_state$max_val)
  K <- model_state$n_components
  prior_obs <- alpha / K

  # Calculate marginal probabilities ignoring NAs
  marginal_prob <- numeric(ncol(X_oh))
  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    cols <- ((j - 1) * model_state$max_val + 1):(j * model_state$max_val)
    if (any(valid)) {
      if (!is.null(weights)) {
        marginal_prob[cols] <- colSums(X_oh[valid, cols, drop=FALSE] * weights[valid]) / sum(weights[valid])
      } else {
        marginal_prob[cols] <- colMeans(X_oh[valid, cols, drop=FALSE])
      }
    }
  }

  pis <- t(resp) %*% X_oh
  pis <- sweep(pis, 2, prior_obs * marginal_prob, "+")

  sum_resp <- colSums(resp) + prior_obs
  pis <- sweep(pis, 1, sum_resp, "/")

  # Re-normalize just to prevent floating point drift
  for (j in seq_len(ncol(X))) {
    cols <- ((j - 1) * model_state$max_val + 1):(j * model_state$max_val)
    pis[, cols] <- sweep(pis[, cols, drop=FALSE], 1, rowSums(pis[, cols, drop=FALSE]), "/")
  }

  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.multinoulli <- function(model_state, X) {
  X_oh <- one_hot(X, model_state$max_val)
  log_eps <- matrix(0, nrow = nrow(X), ncol = model_state$n_components)

  for (c in seq_len(model_state$n_components)) {
    pi_c <- pmax(pmin(model_state$parameters$pis[c, ], 1 - 1e-15), 1e-15)
    log_eps[, c] <- rowSums(sweep(X_oh, 2, log(pi_c), "*"))
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.multinoulli <- function(model_state) {
  n_features <- ncol(model_state$parameters$pis) / model_state$max_val
  return(model_state$n_components * n_features * (model_state$max_val - 1))
}

# ------------------------------------------------------------------------------
# 4. Multinoulli NaN S3 Methods
# ------------------------------------------------------------------------------
#' @exportS3Method init_params multinoulli_nan
init_params.multinoulli_nan <- init_params.multinoulli
#' @exportS3Method n_parameters multinoulli_nan
n_parameters.multinoulli_nan <- n_parameters.multinoulli

# The base multinoulli m-step now cleanly handles NAs via the one_hot zeroing
# and the targeted valid-row marginal calculations!
#' @exportS3Method m_step multinoulli_nan
m_step.multinoulli_nan <- m_step.multinoulli
#' @exportS3Method log_likelihood multinoulli_nan
log_likelihood.multinoulli_nan <- log_likelihood.multinoulli
