# ==============================================================================
# S3 Core EM Algorithm
# ==============================================================================

# Helper to check if structural model contains a covariate
has_covariate <- function(sm) {
  if (is.null(sm)) return(FALSE)
  if (inherits(sm, "covariate")) return(TRUE)
  if (inherits(sm, "nested")) {
    return(any(sapply(sm$models, function(m) inherits(m, "covariate"))))
  }
  return(FALSE)
}

# Helper to initialize random responsibilities
initialize_resp <- function(n_samples, n_components) {
  resp <- matrix(runif(n_samples * n_components), nrow = n_samples)
  resp <- sweep(resp, 1, rowSums(resp), "/")
  return(resp)
}

# E-step: Calculate responsibilities
e_step <- function(model_state, X, Y = NULL) {
  # 1. Measurement model likelihood
  log_prob <- log_likelihood(model_state$mm, X)
  
  # 2. Add Structural model likelihood (if active)
  if (!is.null(Y) && !is.null(model_state$sm)) {
    log_prob <- log_prob + log_likelihood(model_state$sm, Y)
  }
  
  # 3. Add Marginal Prior ONLY if there is no covariate predicting the classes
  if (!has_covariate(model_state$sm)) {
    # Using 1e-15 here just in case, though the Bayesian prior prevents 0s naturally now
    log_weights <- log(model_state$weights + 1e-15) 
    log_prob <- sweep(log_prob, 2, log_weights, "+")
  }
  
  # 4. Log-sum-exp trick for stability
  log_prob_norm <- apply(log_prob, 1, function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  })
  
  log_resp <- sweep(log_prob, 1, log_prob_norm, "-")
  return(list(log_resp = log_resp, log_prob_norm = log_prob_norm))
}

# M-step: Update model parameters
m_step_core <- function(model_state, X, Y, log_resp, alpha = 1.0) {
  resp <- exp(log_resp)
  
  # --- BAYESIAN PRIOR (Class Weights) ---
  K <- model_state$n_components
  prior_obs <- alpha / K
  
  nk <- colSums(resp)
  nk_prior <- nk + prior_obs
  model_state$weights <- nk_prior / sum(nk_prior)
  
  model_state$mm <- m_step(model_state$mm, X, resp)
  
  if (!is.null(Y) && !is.null(model_state$sm)) {
    model_state$sm <- m_step(model_state$sm, Y, resp)
  }
  
  return(model_state)
}

# Run the EM loop for a single random initialization
fit_single_init <- function(model_state, X, Y, max_iter = 1000, abs_tol = 1e-3, rel_tol = 1e-3) {
  n_samples <- nrow(X)
  
  model_state$weights <- rep(1 / model_state$n_components, model_state$n_components)
  
  # Initialize parameters directly
  model_state$mm <- init_params(model_state$mm, X, NULL)
  if (!is.null(Y) && !is.null(model_state$sm)) {
    model_state$sm <- init_params(model_state$sm, Y, NULL)
  }
  
  # Initialize scalar trackers for convergence
  prev_total_ll <- -Inf
  converged <- FALSE
  n_iter <- 0
  
  for (iter in 1:max_iter) {
    # 1. E-STEP
    e_res <- e_step(model_state, X, Y)
    
    # Store the FULL VECTOR in the model state (for weights/BIC)
    log_prob_vector <- e_res$log_prob_norm
    
    # Calculate the SCALAR TOTAL LL for the convergence check
    # This uses the weights we added to support survey data!
    current_total_ll <- sum(model_state$sample_weights * log_prob_vector)
    
    # 2. CONVERGENCE CHECK (Using scalars)
    if (iter > 1) {
      change <- current_total_ll - prev_total_ll
      
      if (!is.na(change) && (abs(change) < abs_tol || abs(change / max(abs(prev_total_ll), 1e-9)) < rel_tol)) {
        converged <- TRUE
        n_iter <- iter
        model_state$lower_bound <- log_prob_vector # Final vector storage
        break
      }
    }
    
    # 3. M-STEP: Update parameters
    model_state <- m_step_core(model_state, X, Y, e_res$log_resp)
    
    # Prepare for next iteration
    prev_total_ll <- current_total_ll
    n_iter <- iter
  }
  
  model_state$converged <- converged
  model_state$n_iter <- n_iter
  model_state$lower_bound <- log_prob_vector # Ensure the vector is returned
  model_state$log_resp <- e_res$log_resp
  
  return(model_state)
}

# Multi-start EM
fit_em <- function(model_state, X, Y, n_init = 1, max_iter = 1000, random_state = NULL) {
  best_model <- NULL
  # We use a scalar tracker for the comparison
  best_total_ll <- -Inf 
  
  for (init in 1:n_init) {
    if (!is.null(random_state)) set.seed(random_state + init)
    
    # Run the single initialization
    fitted_state <- fit_single_init(model_state, X, Y, max_iter = max_iter)
    
    # 1. CALCULATE SCALAR SUM for comparison
    # fitted_state$lower_bound is the vector
    # current_ll is the single scalar sum
    current_ll <- sum(fitted_state$sample_weights * fitted_state$lower_bound)
    
    # 2. COMPARE SCALARS
    # This prevents the 'length = 3000' error
    if (is.null(best_model) || current_ll > best_total_ll) {
      best_total_ll <- current_ll
      best_model <- fitted_state
    }
  }
  
  return(best_model)
}