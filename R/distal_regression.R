# ==============================================================================
# S3 Distal Regression Model (Simultaneous Y ~ X + Z, class-varying slopes)
# ==============================================================================
# distal_one_hot() and distal_forward() are defined in utils.R

distal_regression_model <- function(n_components, tol = 1e-4, max_iter = 100,
                                    method = "newton-raphson") {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                method = method, parameters = list())
  class(state) <- c("distal_regression", "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# Internal helper: validate and normalise the categorical outcome column.
# ------------------------------------------------------------------------------
.validate_distal_Y <- function(Y_col, context = "distal_regression") {
  Y <- as.numeric(Y_col)

  # Reject negative values (Bug 5 fix)
  if (any(!is.na(Y) & Y < 0))
    stop(sprintf(
      paste0("%s: Y outcome column contains negative values (%s). ",
             "Outcomes must be 0-indexed (0, 1, ...) or 1-indexed (1, 2, ...)."),
      context,
      paste(sort(unique(Y[!is.na(Y) & Y < 0]))[1:min(3, sum(!is.na(Y) & Y < 0))],
            collapse = ", ")
    ))

  # Shift 0-indexed to 1-indexed
  if (!all(is.na(Y)) && min(Y, na.rm = TRUE) == 0) Y <- Y + 1

  # Guard against degenerate constant outcome (Bug 4 fix)
  M <- max(Y, na.rm = TRUE)
  if (M <= 1) {
    warning(sprintf(
      paste0("%s: outcome Y has only one unique category (M = %d). ",
             "A regression model cannot be estimated. ",
             "Returning the model state unchanged."),
      context, M
    ))
    return(NULL)  # signals callers to skip the NR loop
  }

  return(Y)
}

#' @exportS3Method
init_params.distal_regression <- function(model_state, X, resp,
                                          random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)

  Y <- .validate_distal_Y(X[, 1], "distal_regression init_params")
  if (is.null(Y)) {
    # Degenerate constant outcome: store a placeholder so downstream code
    # doesn't crash on a missing max_val.
    model_state$max_val <- 1L
    model_state$parameters$betas <- array(0, dim = c(model_state$n_components, 0, 0))
    return(model_state)
  }

  Z <- cbind(1, X[, -1, drop = FALSE])

  K <- model_state$n_components
  D <- ncol(Z)
  M <- max(Y, na.rm = TRUE)

  model_state$max_val <- M
  model_state$parameters$betas <-
    array(rnorm(K * (M - 1) * D, 0, 1.0), dim = c(K, M - 1, D))
  return(model_state)
}

#' @exportS3Method
m_step.distal_regression <- function(model_state, X, resp, weights = NULL, ...) {
  Y <- .validate_distal_Y(X[, 1], "distal_regression m_step")
  if (is.null(Y)) return(model_state)  # constant outcome — skip (Bug 4 fix)

  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z     <- cbind(1, Z_raw)

  M       <- model_state$max_val
  K       <- model_state$n_components
  D       <- ncol(Z)
  Y_oh    <- distal_one_hot(Y, M)
  valid_Y <- !is.na(Y)

  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")
  model_state$parameters$hessians <- list()

  for (k in seq_len(K)) {
    W_k        <- resp[valid_Y, k]
    Z_valid    <- Z[valid_Y, , drop = FALSE]
    Y_oh_valid <- Y_oh[valid_Y, , drop = FALSE]
    Nk_eff     <- sum(W_k)

    beta_k <- model_state$parameters$betas[k, , ]
    # t(as.matrix(vector)) transposes incorrectly when D=1: gives 1x(M-1)
    # instead of (M-1)x1. Use matrix() with explicit dims always.
    if (!is.matrix(beta_k))
      beta_k <- matrix(beta_k, nrow = M - 1, ncol = D)

    # ----------------------------------------------------------------
    # Guard: skip optimization for near-empty classes.
    # A class with negligible effective weight has no information to
    # estimate its regression parameters; attempting NR on such a class
    # produces a near-singular Hessian and divergent parameter updates.
    # ----------------------------------------------------------------
    if (Nk_eff > 1.0) {
      for (iter in seq_len(model_state$max_iter)) {
        P <- distal_forward(Z_valid, beta_k)

        # Gradient
        G <- matrix(0, nrow = M - 1, ncol = D)
        for (m in 2:M) {
          diff       <- W_k * (Y_oh_valid[, m] - P[, m])
          G[m - 1, ] <- t(Z_valid) %*% diff
        }
        G_flat <- as.vector(t(G))
        if (max(abs(G_flat)) < model_state$tol) break

        # Hessian for NR step (with ridge stabilisation)
        H_nr <- matrix(0, nrow = (M - 1) * D, ncol = (M - 1) * D)
        for (m in 2:M) {
          for (m_prime in 2:M) {
            kronecker <- if (m == m_prime) 1 else 0
            r     <- W_k * P[, m] * (kronecker - P[, m_prime])
            sub_H <- -t(Z_valid) %*% sweep(Z_valid, 1, r, "*")
            idx_m       <- ((m - 2) * D + 1):((m - 1) * D)
            idx_m_prime <- ((m_prime - 2) * D + 1):((m_prime - 1) * D)
            H_nr[idx_m, idx_m_prime] <- sub_H
          }
        }
        diag(H_nr) <- diag(H_nr) - 1e-4

        # Step with clipping — prevents divergence for small classes
        update_step <- pinv(-H_nr) %*% G_flat
        max_step    <- 1.0
        if (max(abs(update_step)) > max_step)
          update_step <- update_step * (max_step / max(abs(update_step)))

        beta_flat <- as.vector(t(beta_k)) + update_step
        beta_k    <- t(matrix(beta_flat, nrow = D, ncol = M - 1))
      }
    }

    # ----------------------------------------------------------------
    # Always recompute H at the FINAL converged beta_k.
    # Bug fix: previously H was only assigned inside the loop, after the
    # gradient check — meaning it was never set when the loop broke early
    # (convergence at first iteration).  The stale value from the previous
    # class's loop was silently used instead, producing identical Hessians
    # across classes.  Computing H here, outside the loop and without the
    # ridge penalty, gives the true observed Fisher information at the MLE.
    # ----------------------------------------------------------------
    P_final  <- distal_forward(Z_valid, beta_k)
    H_final  <- matrix(0, nrow = (M - 1) * D, ncol = (M - 1) * D)
    for (m in 2:M) {
      for (m_prime in 2:M) {
        kronecker <- if (m == m_prime) 1 else 0
        r     <- W_k * P_final[, m] * (kronecker - P_final[, m_prime])
        sub_H <- -t(Z_valid) %*% sweep(Z_valid, 1, r, "*")
        idx_m       <- ((m - 2) * D + 1):((m - 1) * D)
        idx_m_prime <- ((m_prime - 2) * D + 1):((m_prime - 1) * D)
        H_final[idx_m, idx_m_prime] <- sub_H
      }
    }

    model_state$parameters$betas[k, , ]  <- beta_k
    model_state$parameters$hessians[[k]] <- H_final
  }
  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_regression <- function(model_state, X, ...) {
  Y <- .validate_distal_Y(X[, 1], "distal_regression log_likelihood")
  if (is.null(Y)) {
    # Constant outcome: contribute zero log-likelihood (uninformative)
    return(matrix(0, nrow = nrow(X), ncol = model_state$n_components))
  }

  Z_raw <- impute_covariates(X[, -1, drop = FALSE])
  Z     <- cbind(1, Z_raw)

  K     <- model_state$n_components
  M     <- model_state$max_val
  D     <- dim(model_state$parameters$betas)[3]  # must be defined before the loop
  n     <- nrow(X)
  valid <- !is.na(Y)

  log_eps <- matrix(0, nrow = n, ncol = K)

  for (k in seq_len(K)) {
    beta_k <- model_state$parameters$betas[k, , ]
    # t(as.matrix(vector)) transposes incorrectly when D=1: gives 1x(M-1)
    # instead of (M-1)x1. Use matrix() with explicit dims always.
    if (!is.matrix(beta_k))
      beta_k <- matrix(beta_k, nrow = M - 1, ncol = D)

    P <- pmax(pmin(distal_forward(Z, beta_k), 1 - 1e-15), 1e-15)

    prob_obs <- numeric(n)
    if (any(valid)) prob_obs[valid] <- P[cbind(which(valid), Y[valid])]
    prob_obs[!valid] <- 1

    log_eps[, k] <- log(prob_obs)
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.distal_regression <- function(model_state, ...) {
  K <- model_state$n_components
  D <- dim(model_state$parameters$betas)[3]
  M <- model_state$max_val
  return(K * (M - 1) * D)
}
