# ==============================================================================
# S3 Nested Model Dispatcher
# ==============================================================================

#' @exportS3Method
init_params.nested <- function(model_state, X, resp, random_state = NULL) {
  start_col <- 1
  for (name in names(model_state$models)) {
    n_cols <- model_state$columns_per_model[name]
    end_col <- start_col + n_cols - 1

    X_sub <- X[, start_col:end_col, drop = FALSE]
    model_state$models[[name]] <- init_params(model_state$models[[name]], X_sub, resp, random_state)

    start_col <- end_col + 1
  }
  return(model_state)
}

#' @exportS3Method
m_step.nested <- function(model_state, X, resp, ...) {
  start_col <- 1
  for (name in names(model_state$models)) {
    n_cols <- model_state$columns_per_model[name]
    end_col <- start_col + n_cols - 1

    X_sub <- X[, start_col:end_col, drop = FALSE]
    model_state$models[[name]] <- m_step(model_state$models[[name]], X_sub, resp, ...)

    start_col <- end_col + 1
  }
  return(model_state)
}

#' @exportS3Method
log_likelihood.nested <- function(model_state, X) {
  n_samples <- nrow(X)
  log_eps <- matrix(0, nrow = n_samples, ncol = model_state$n_components)

  start_col <- 1
  for (name in names(model_state$models)) {
    n_cols <- model_state$columns_per_model[name]
    end_col <- start_col + n_cols - 1

    X_sub <- X[, start_col:end_col, drop = FALSE]
    log_eps <- log_eps + log_likelihood(model_state$models[[name]], X_sub)

    start_col <- end_col + 1
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.nested <- function(model_state) {
  n <- 0
  for (name in names(model_state$models)) {
    n <- n + n_parameters(model_state$models[[name]])
  }
  return(n)
}
