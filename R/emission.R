# ==============================================================================
# S3 Generics for Emission Models (Replacing R6 Base Class)
# ==============================================================================

# Initialize parameters for an emission model
init_params <- function(model_state, X, resp, ...) {
  UseMethod("init_params")
}

# Perform the M-step (update parameters based on responsibilities)
m_step <- function(model_state, X, resp, ...) {
  UseMethod("m_step")
}

# Calculate the log-likelihood (E-step)
log_likelihood <- function(model_state, X, ...) {
  UseMethod("log_likelihood")
}

# Count the number of free parameters
n_parameters <- function(model_state, ...) {
  UseMethod("n_parameters")
}
