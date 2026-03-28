# ==============================================================================
# S3 Covariate Model - FINAL (Firth Estimation + Unpenalized Inference)
# ==============================================================================

covariate_model <- function(n_components, tol = 1e-6, max_iter = 500, intercept = TRUE) {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                intercept = intercept, parameters = list())
  class(state) <- c("covariate", "emission")
  return(state)
}

#' @exportS3Method
m_step.covariate <- function(model_state, X, resp, weights = NULL, ...) {
  X_mat <- if(model_state$intercept) cbind(1, as.matrix(X)) else as.matrix(X)
  K <- model_state$n_components
  D <- ncol(X_mat)
  w_vec <- if(!is.null(weights)) weights else rep(1, nrow(X_mat))

  # ============================================================================
  # 1. ESTIMATION PASS (Data Augmentation)
  # ============================================================================
  # We add a tiny "ghost" observation to every class at the mean of X.
  # This makes Complete Separation mathematically impossible during optimization.
  pseudo_X <- matrix(colMeans(X_mat), nrow = K, ncol = D, byrow = TRUE)
  pseudo_resp <- diag(K)
  pseudo_w <- rep(0.01, K)

  X_aug <- rbind(X_mat, pseudo_X)
  resp_aug <- rbind(resp, pseudo_resp)
  w_aug <- c(w_vec, pseudo_w)

  nll_func <- function(pars) {
    B <- rbind(matrix(pars, K-1, D, byrow=TRUE), 0)
    logits <- X_aug %*% t(B)
    logits <- pmax(pmin(logits, 50), -50)
    max_l <- apply(logits, 1, max)
    prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))
    -sum(w_aug * rowSums(resp_aug * log(prob + 1e-15)))
  }

  # Analytical gradient of nll_func wrt pars.
  # For multinomial logistic regression:
  #   ∂nll/∂B[k,d] = Σ_i w_i (prob[i,k] - resp[i,k]) * X[i,d]  for k < K
  # Vectorised: grad[k, :] = t(X_aug) %*% (w_aug * (prob[,k] - resp_aug[,k]))
  nll_grad <- function(pars) {
    B <- rbind(matrix(pars, K-1, D, byrow=TRUE), 0)
    logits <- X_aug %*% t(B)
    logits <- pmax(pmin(logits, 50), -50)
    max_l <- apply(logits, 1, max)
    prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))
    # residual: (prob - resp_aug) weighted by w_aug, for free classes k = 1..K-1
    resid <- sweep(prob[, seq_len(K-1), drop=FALSE] -
                     resp_aug[, seq_len(K-1), drop=FALSE], 1, w_aug, "*")
    as.vector(t(X_aug) %*% resid)   # D × (K-1), flattened row-major
  }

  fit <- optim(par = rep(0, (K-1)*D), fn = nll_func, gr = nll_grad,
               method = "BFGS")
  beta_final <- rbind(matrix(fit$par, K-1, D, byrow=TRUE), 0)

  # ============================================================================
  # 2. INFERENCE PASS (Original Data Only)
  # ============================================================================
  # We throw away the pseudo-data and calculate the Hessian strictly on the
  # original X_mat and w_vec to ensure Standard Errors are not artificially shrunk.
  logits <- X_mat %*% t(beta_final)
  logits <- pmax(pmin(logits, 50), -50)
  max_l <- apply(logits, 1, max)
  prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))

  H <- matrix(0, (K-1)*D, (K-1)*D)

  for(k in seq_len(K-1)) {
    for(j in seq_len(K-1)) {
      W <- prob[,k] * ((if(k==j) 1 else 0) - prob[,j]) * w_vec
      H_kj <- -t(X_mat) %*% sweep(X_mat, 1, W, "*")

      H[((k-1)*D+1):(k*D), ((j-1)*D+1):(j*D)] <- H_kj
    }
  }

  # Pad with zeros for the anchor class (identifiability constraint)
  H_full <- matrix(0, K*D, K*D)
  if ((K-1)*D > 0)
    H_full[1:((K-1)*D), 1:((K-1)*D)] <- H
  diag(H_full)[((K-1)*D+1):(K*D)] <- -1e8

  model_state$parameters$beta <- beta_final
  model_state$parameters$hessian <- H_full
  return(model_state)
}

#' @exportS3Method
init_params.covariate <- function(model_state, X, resp, ...) {
  D <- ncol(X) + as.integer(model_state$intercept)
  model_state$parameters$beta <- matrix(0, model_state$n_components, D)
  return(model_state)
}

#' @exportS3Method
log_likelihood.covariate <- function(model_state, X, ...) {
  # 1. Prepare Matrix
  X_mat <- if(model_state$intercept) cbind(1, as.matrix(X)) else as.matrix(X)

  # 2. Compute Raw Logits
  logits <- X_mat %*% t(model_state$parameters$beta)

  # 3. Use the Log-Sum-Exp trick for the denominator
  # This uses your existing helper in utils.R
  log_denominator <- logsumexp(logits, MARGIN = 1)

  # 4. Compute Log-Probabilities: Log(Prob) = Logits - Log(Sum_Exp_Logits)
  # sweep subtracts the log_denominator from each row of logits
  # sweep() can silently drop the dim attribute in ALTREP system,
  # causing non-conformable-arrays when log_prob is added to the mm log-likelihood
  # in e_step. matrix() re-attaches explicit dimensions unconditionally.
  log_prob <- matrix(
    sweep(logits, 1, log_denominator, "-"),
    nrow = nrow(logits), ncol = ncol(logits)
  )

  return(log_prob)
}

#' @exportS3Method
n_parameters.covariate <- function(model_state, ...) {
  return((nrow(model_state$parameters$beta) - 1) * ncol(model_state$parameters$beta))
}
