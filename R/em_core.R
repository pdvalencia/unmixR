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

  # 3. Add Marginal Prior ONLY if the covariate model is not actively providing
  # class probabilities.  A covariate SM is "active" only when Y is non-NULL AND
  # the SM is a covariate type — in that case log_likelihood(sm, Y) already
  # encodes P(class | z_i) so adding log_weights would double-count.
  # When Y = NULL (e.g. during step-1 of 3-step estimation), the covariate SM
  # contributes nothing and marginal class weights MUST still be applied.
  covariate_active <- !is.null(Y) && has_covariate(model_state$sm)
  if (!covariate_active) {
    log_weights <- log(model_state$weights + 1e-15)
    log_prob    <- sweep(log_prob, 2, log_weights, "+")
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

# L-BFGS refinement after EM convergence — Penalised Maximum Likelihood (PM).
#
# EM converges to a Q-function fixed point, not necessarily the PM optimum.
# L-BFGS steps on the full penalised log-posterior climb past that fixed point.
#
# The objective matches LG's PM formulation (technical guide §7.3–7.4):
#   log P(ϑ) = log L(X; ϑ) + log p(ϑ)
# where the Dirichlet priors (Bayes constants α1 = α2 = 1, LG defaults) are:
#   Weights  : (α1/K)   · Σ_k  log π_k
#   Bernoulli: (α2/K)   · Σ_k Σ_j  [ π̂_j · log π_kj + (1-π̂_j) · log(1-π_kj) ]
# with π̂_j the observed marginal probability of item j (LG's "conservative null
# model"). Gaussian models carry no Dirichlet prior in the EM M-step (gaussian_diag
# uses a hard-floor regularisation instead), so only the weight prior is added.
#
# Using the prior instead of box constraints is the correct approach: it lets
# the gradient pull parameters freely while the prior provides soft penalisation
# proportional to the evidence — exactly what LG does in its NR phase (eq. 20–21).
#
# Supported: bernoulli, bernoulli_nan, gaussian_diag, gaussian_unit.
# No-op for nested models and other types.
#
# Parameterisation (unconstrained):
#   bernoulli     : logit(pis) [K×J]  + log-ratio weights [K-1]
#   gaussian_diag : means [K×J] + log(sd) [K×J] + log-ratio weights [K-1]
#   gaussian_unit : means [K×J] + log-ratio weights [K-1]
#
refine_lbfgs <- function(model_state, X, Y = NULL, max_iter = 500) {
  mm_type <- class(model_state$mm)[1]
  supported <- c("bernoulli", "bernoulli_nan", "gaussian_diag", "gaussian_unit")
  if (!mm_type %in% supported) return(model_state)
  if (inherits(model_state$mm, "nested")) return(model_state)
  # K=1 has no weight parameters; the M-step already gives the exact analytic
  # solution (item marginals), so L-BFGS is a no-op and the K-2 index arithmetic
  # below produces an out-of-bounds sequence that triggers a sweep() warning.
  if (model_state$n_components == 1L) return(model_state)

  K  <- model_state$n_components
  J  <- ncol(X)
  sw <- model_state$sample_weights

  # Observed marginal probabilities — the "conservative null model" LG uses
  # as the centre of the Dirichlet prior for binary items (guide §7.3).
  marginal <- colMeans(X, na.rm = TRUE)
  marginal <- pmax(pmin(marginal, 1 - 1e-7), 1e-7)

  # ── Pack initial parameters ────────────────────────────────────────────────
  wts <- pmax(model_state$weights, 1e-15)
  log_ratio_w <- log(wts[-K] / wts[K])   # K-1 free weight params (last anchored = 0)

  if (mm_type %in% c("bernoulli", "bernoulli_nan")) {
    pis  <- pmax(pmin(model_state$mm$parameters$pis, 1 - 1e-7), 1e-7)
    par0 <- c(qlogis(as.vector(pis)), log_ratio_w)   # K*J + K-1

  } else if (mm_type == "gaussian_diag") {
    means <- as.vector(model_state$mm$parameters$means)
    sds   <- sqrt(as.vector(model_state$mm$parameters$covariances))
    par0  <- c(means, log(pmax(sds, 1e-7)), log_ratio_w)   # 2*K*J + K-1

  } else {  # gaussian_unit
    par0 <- c(as.vector(model_state$mm$parameters$means), log_ratio_w)   # K*J + K-1
  }

  n_obs <- sum(sw)   # effective sample size (supports survey weights)

  # ── Penalised log-posterior + analytical gradient ──────────────────────────
  # Providing an analytical gradient to optim() eliminates the O(p) function
  # evaluations per iteration that numerical finite-difference gradient
  # estimation requires (p = K*J + K-1 ≈ 47 for K=4, J=11).  Without it,
  # L-BFGS spends ~50 * (p+1) ≈ 5700 function calls per fit; with it, the
  # inner loop reduces to ~50 * 4 ≈ 200 calls — a ~25x speedup that makes
  # refine_lbfgs fast at any sample size.
  #
  # The gradient derivations follow from the chain rule on the PM objective
  # PM(θ) = log L(X;θ) + log p(θ) and the reparameterisations:
  #   · Bernoulli pis on logit scale:  ∂PM/∂logit(π_kj) = Σ_i sw_i R_ik (x_ij − π_kj) + (m̂_j − π_kj)/K
  #   · Class weights on log-ratio:   ∂PM/∂lr_k = Σ_i sw_i (R_ik − w_k) + 1/K − w_k  (k < K)
  #   · Gaussian means (unit var):     ∂PM/∂μ_kj = Σ_i sw_i R_ik (x_ij − μ_kj)
  #   · Gaussian means (diag):         ∂PM/∂μ_kj = Σ_i sw_i R_ik (x_ij − μ_kj) / σ_kj²
  #   · Gaussian log-sd (diag):        ∂PM/∂log(σ_kj) = Σ_i sw_i R_ik [(x_ij−μ_kj)²/σ_kj² − 1]

  # Helper: given current params, compute neg-PM value AND gradient in one pass.
  neg_pm_and_grad <- function(par) {

    # ── Decode class weights ─────────────────────────────────────────────────
    n_w  <- K - 1L
    lr   <- c(par[(length(par) - n_w + 1):length(par)], 0)
    lr   <- lr - (max(lr) + log(sum(exp(lr - max(lr)))))
    log_w <- lr
    w_vec <- exp(log_w)

    # ── Decode measurement model ─────────────────────────────────────────────
    if (mm_type %in% c("bernoulli", "bernoulli_nan")) {
      pis_p <- pmax(pmin(matrix(plogis(par[seq_len(K * J)]), K, J), 1-1e-15), 1e-15)
      # Vectorised log-likelihood: n×K  (matrix multiply — no R for-loop over k)
      log_lik <- X %*% t(log(pis_p)) + (1 - X) %*% t(log(1 - pis_p))

    } else if (mm_type == "gaussian_diag") {
      means_p <- matrix(par[seq_len(K * J)],           K, J)
      sds_p   <- matrix(exp(par[(K*J+1):(2*K*J)]),     K, J)
      sds_p   <- pmax(sds_p, 1e-7)
      log_lik <- matrix(0, nrow(X), K)
      for (k in seq_len(K))
        log_lik[, k] <- rowSums(
          dnorm(X, mean = matrix(means_p[k,], nrow(X), J, byrow = TRUE),
                sd   = matrix(sds_p[k,],  nrow(X), J, byrow = TRUE), log = TRUE),
          na.rm = TRUE)

    } else {  # gaussian_unit
      means_p <- matrix(par[seq_len(K * J)], K, J)
      log_lik <- matrix(0, nrow(X), K)
      for (k in seq_len(K))
        log_lik[, k] <- rowSums(
          dnorm(X, mean = matrix(means_p[k,], nrow(X), J, byrow = TRUE),
                sd = 1, log = TRUE), na.rm = TRUE)
    }

    # ── Posterior responsibilities (needed for both value and gradient) ───────
    log_joint <- sweep(log_lik, 2, log_w, "+")
    mx        <- apply(log_joint, 1, max)
    log_norm  <- mx + log(rowSums(exp(sweep(log_joint, 1, mx, "-"))))
    R         <- exp(sweep(log_joint, 1, log_norm, "-"))   # n × K posteriors

    # ── Observed-data log-likelihood (normalised by n_obs for scale-invariance) ─
    obs_ll <- sum(sw * log_norm) / n_obs

    # ── Priors (also normalised by n_obs) ─────────────────────────────────────
    log_prior_w <- (1 / K) * sum(log_w) / n_obs
    log_prior_pis <- if (mm_type %in% c("bernoulli", "bernoulli_nan"))
      sum((marginal / K) %*% t(log(pis_p)) + ((1 - marginal) / K) %*% t(log(1 - pis_p))) / n_obs
    else 0

    val <- -(obs_ll + log_prior_w + log_prior_pis)

    # ── Analytical gradient ───────────────────────────────────────────────────
    swR <- sweep(R, 1, sw, "*")          # sw_i * R_ik,  n × K
    nk  <- colSums(swR)                  # Σ_i sw_i R_ik, length K

    grad <- numeric(length(par))

    if (mm_type %in% c("bernoulli", "bernoulli_nan")) {
      g_pis <- t(swR) %*% X -
        sweep(pis_p, 1, nk, "*") +
        (matrix(marginal, K, J, byrow = TRUE) - pis_p) / K
      grad[seq_len(K * J)] <- as.vector(-g_pis) / n_obs

    } else if (mm_type == "gaussian_diag") {
      g_mu  <- matrix(0, K, J)
      g_lsd <- matrix(0, K, J)
      for (k in seq_len(K)) {
        res_k <- sweep(X, 2, means_p[k,], "-")
        g_mu[k,]  <- colSums(swR[,k] * sweep(res_k, 2, sds_p[k,]^(-2), "*"))
        g_lsd[k,] <- colSums(swR[,k] * sweep(res_k^2, 2, sds_p[k,]^(-2), "*") - nk[k])
      }
      grad[seq_len(K * J)]  <- as.vector(-g_mu)  / n_obs
      grad[(K*J+1):(2*K*J)] <- as.vector(-g_lsd) / n_obs

    } else {  # gaussian_unit
      g_mu <- matrix(0, K, J)
      for (k in seq_len(K))
        g_mu[k,] <- colSums(swR[,k] * sweep(X, 2, means_p[k,], "-"))
      grad[seq_len(K * J)] <- as.vector(-g_mu) / n_obs
    }

    g_w <- nk[seq_len(n_w)] - n_obs * w_vec[seq_len(n_w)] +
      1/K - w_vec[seq_len(n_w)]
    grad[(length(par) - n_w + 1):length(par)] <- -g_w / n_obs

    list(val = val, grad = grad)
  }

  # Wrapper returning just value (for optim fn=)
  neg_pm_val  <- function(par) neg_pm_and_grad(par)$val
  # Wrapper returning just gradient (for optim gr=)
  neg_pm_grad <- function(par) neg_pm_and_grad(par)$grad

  # ── Run L-BFGS with analytical gradient (unconstrained) ────────────────────
  fit <- tryCatch(
    optim(par0, neg_pm_val, gr = neg_pm_grad, method = "L-BFGS-B",
          control = list(maxit = max_iter, factr = 1e7)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence > 1) return(model_state)

  # ── Unpack ────────────────────────────────────────────────────────────────
  par <- fit$par
  n_w <- K - 1L

  lr   <- c(par[(length(par) - n_w + 1):length(par)], 0)
  lr   <- lr - (max(lr) + log(sum(exp(lr - max(lr)))))
  model_state$weights <- exp(lr)

  if (mm_type %in% c("bernoulli", "bernoulli_nan")) {
    model_state$mm$parameters$pis <-
      pmax(pmin(matrix(plogis(par[1:(K*J)]), nrow=K, ncol=J), 1-1e-7), 1e-7)

  } else if (mm_type == "gaussian_diag") {
    model_state$mm$parameters$means       <- matrix(par[1:(K*J)],          nrow=K, ncol=J)
    model_state$mm$parameters$covariances <- matrix(exp(par[(K*J+1):(2*K*J)]), nrow=K, ncol=J)^2

  } else {  # gaussian_unit
    model_state$mm$parameters$means <- matrix(par[1:(K*J)], nrow=K, ncol=J)
  }

  # Re-run E-step with the refined parameters so log_resp and lower_bound
  # reflect the true post-refinement posteriors.
  e_res <- e_step(model_state, X, Y)
  model_state$log_resp    <- e_res$log_resp
  model_state$lower_bound <- e_res$log_prob_norm

  return(model_state)
}

# Run the EM loop for a single random initialization
fit_single_init <- function(model_state, X, Y, max_iter = 1000,
                            abs_tol = 1e-3, rel_tol = 1e-3, refine = TRUE) {
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

  # L-BFGS refinement per restart so all restarts compete on PM likelihood.
  # With analytical gradients this is fast (~0.3s at n=5000), so the n_init-fold
  # cost is acceptable and selection bias from EM-only ranking is avoided.
  # Skipped when refine = FALSE (BLRT bootstrap replicates).
  if (isTRUE(refine)) model_state <- refine_lbfgs(model_state, X, Y)

  return(model_state)
}

# Multi-start EM
fit_em <- function(model_state, X, Y, n_init = 1, max_iter = 1000,
                   random_state = NULL, refine = TRUE) {
  best_model <- NULL
  best_total_ll <- -Inf

  for (init in 1:n_init) {
    if (!is.null(random_state)) set.seed(random_state + init)

    # Each restart includes L-BFGS refinement so selection is by PM likelihood.
    fitted_state <- fit_single_init(model_state, X, Y, max_iter = max_iter,
                                    refine = refine)

    current_ll <- sum(fitted_state$sample_weights * fitted_state$lower_bound)
    if (is.null(best_model) || current_ll > best_total_ll) {
      best_total_ll <- current_ll
      best_model    <- fitted_state
    }
  }

  return(best_model)
}
