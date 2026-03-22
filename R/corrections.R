# ==============================================================================
# S3 3-Step Corrections (BCH and ML)
# ==============================================================================

get_modal_resp <- function(resp) {
  n <- nrow(resp)
  K <- ncol(resp)
  modal <- matrix(0, nrow = n, ncol = K)
  assigned_classes <- max.col(resp, ties.method = "first")
  for (i in 1:n) {
    modal[i, assigned_classes[i]] <- 1
  }
  return(modal)
}

# Apply the BCH Correction
fit_bch <- function(model_state, X, Y) {

  # Warn if BCH is used with outcome types for which ML is preferred
  # (Vermunt, 2010; Bakk, Tekle & Vermunt, 2013)
  if (inherits(model_state$sm, c("covariate", "distal_regression", "distal_pooled"))) {
    warning(paste(
      "BCH correction is not recommended for covariates or categorical distal outcomes.",
      "Vermunt (2010) and Bakk et al. (2013) show ML correction is preferred in this case.",
      "Consider correction = 'ML' instead."
    ))
  }

  weights <- model_state$sample_weights
  e_res <- e_step(model_state, X, NULL)
  resp <- exp(e_res$log_resp)

  modal_resp <- get_modal_resp(resp)

  # C[j, k] = P(Assigned=j | True=k)
  C <- t(modal_resp) %*% (resp * weights)
  C <- sweep(C, 2, colSums(resp * weights), "/")

  # According to Bakk et al. (2013), the weights are the columns of the inverse.
  D <- t(pinv(C))

  bch_resp <- modal_resp %*% D

  # Clip negative weights — BCH can produce negatives for poorly-separated
  # classes. Truncating to 0 (rather than abs()) is conservative but honest:
  # it discards ambiguous observations rather than misrepresenting them.
  bch_resp <- pmax(bch_resp, 0)

  model_state$sm <- init_params(model_state$sm, Y, bch_resp)
  model_state$sm <- m_step(model_state$sm, Y, bch_resp)

  return(model_state)
}

# Apply the Maximum Likelihood (ML) Correction (Robust to Missing Covariates)
fit_ml <- function(model_state, X, Y, max_iter = 1000, abs_tol = 1e-3, rel_tol = 1e-3) {

  # Warn if ML is used with continuous outcomes, for which BCH is preferred
  # (Bakk & Vermunt, 2016)
  if (inherits(model_state$sm, c("distal_continuous", "distal_continuous_regression"))) {
    warning(paste(
      "ML correction is not recommended for continuous distal outcomes.",
      "Bakk & Vermunt (2016) show BCH is preferred in this case,",
      "as ML requires strong distributional assumptions.",
      "Consider correction = 'BCH' instead."
    ))
  }

  keep <- complete.cases(Y)
  X_clean <- X[keep, , drop = FALSE]
  Y_clean <- Y[keep, , drop = FALSE]
  w_clean <- model_state$sample_weights[keep]

  if (nrow(X_clean) == 0) stop("All rows have missing covariates.")

  e_res_clean <- e_step(model_state, X_clean, NULL)
  resp_clean <- exp(e_res_clean$log_resp)

  model_state$sm <- init_params(model_state$sm, Y_clean, resp_clean)
  model_state$sm <- m_step(model_state$sm, Y_clean, resp_clean)

  prev_total_ll <- -Inf

  for (iter in 1:max_iter) {
    e_res <- e_step(model_state, X_clean, Y_clean)

    current_total_ll <- sum(w_clean * e_res$log_prob_norm)

    if (iter > 1) {
      change <- current_total_ll - prev_total_ll
      denom <- max(abs(prev_total_ll), 1e-9)
      if (is.na(change) || abs(change) < abs_tol || abs(change / denom) < rel_tol) {
        break
      }
    }

    joint_resp <- exp(e_res$log_resp)
    nk <- colSums(joint_resp * w_clean)
    model_state$weights <- nk / sum(nk)

    model_state$sm <- m_step(model_state$sm, Y_clean, joint_resp, weights = w_clean)
    prev_total_ll <- current_total_ll
  }

  # Final re-alignment across the full dataset.
  # e_step returns log_resp and log_prob_norm but not log_prob directly.
  # Reconstruct: log_prob = log_resp + log_prob_norm (row-wise), since
  # log_resp = log_prob - log_prob_norm by definition.
  e_res_full <- e_step(model_state, X, NULL)
  log_prob_full_mm <- sweep(e_res_full$log_resp, 1, e_res_full$log_prob_norm, "+")

  K <- model_state$n_components
  sm_ll_full <- matrix(0, nrow = nrow(X), ncol = K)
  if (any(keep)) {
    sm_ll_full[keep, ] <- log_likelihood(model_state$sm, Y_clean)
  }

  log_prob_full <- log_prob_full_mm + sm_ll_full
  max_log_prob_full <- apply(log_prob_full, 1, max)
  log_resp_full <- sweep(log_prob_full, 1, max_log_prob_full, "-")
  log_resp_sums <- log(rowSums(exp(log_resp_full)))

  model_state$log_resp <- sweep(log_resp_full, 1, log_resp_sums, "-")
  model_state$lower_bound <- max_log_prob_full + log_resp_sums

  return(model_state)
}
