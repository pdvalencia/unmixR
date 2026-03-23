# ==============================================================================
# S3 User Wrappers and Pipeline Tools (LCA/LPA Mixture Engine)
# ==============================================================================

sort_model_classes <- function(model_state) {
  K <- model_state$n_components
  if (K <= 1) return(model_state)

  new_order <- order(model_state$weights, decreasing = TRUE)
  model_state$weights <- model_state$weights[new_order]

  # --- Sort flat measurement model parameters ---
  if (!is.null(model_state$mm$parameters[["pis"]])) {
    model_state$mm$parameters$pis <-
      model_state$mm$parameters$pis[new_order, , drop = FALSE]
  }
  # Gaussian LPA: also sort means and covariances
  if (!is.null(model_state$mm$parameters[["means"]])) {
    model_state$mm$parameters$means <-
      model_state$mm$parameters$means[new_order, , drop = FALSE]
  }
  if (!is.null(model_state$mm$parameters[["covariances"]])) {
    model_state$mm$parameters$covariances <-
      model_state$mm$parameters$covariances[new_order, , drop = FALSE]
  }

  # --- Sort nested measurement model sub-model parameters ---
  # The flat-parameter block above only touches model_state$mm$parameters, which
  # is empty for nested models.  Sub-model parameters live one level deeper at
  # model_state$mm$models[[name]]$parameters and must be sorted independently.
  if (inherits(model_state$mm, "nested")) {
    for (name in names(model_state$mm$models)) {
      sub <- model_state$mm$models[[name]]
      if (!is.null(sub$parameters[["pis"]]))
        sub$parameters$pis <- sub$parameters$pis[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["means"]]))
        sub$parameters$means <- sub$parameters$means[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["covariances"]]))
        sub$parameters$covariances <-
          sub$parameters$covariances[new_order, , drop = FALSE]
      model_state$mm$models[[name]] <- sub
    }
  }

  # --- Sort structural model parameters ---
  if (!is.null(model_state$sm)) {

    sort_sm_params <- function(sm) {
      if (!is.null(sm$parameters[["beta"]])) {
        sm$parameters$beta <- sm$parameters$beta[new_order, , drop = FALSE]
        if (!is.null(sm$parameters[["hessian"]])) {
          H       <- sm$parameters$hessian
          D       <- ncol(sm$parameters$beta)
          idx_map <- as.vector(sapply(new_order, function(k) ((k-1)*D + 1):(k*D)))
          sm$parameters$hessian <- H[idx_map, idx_map, drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["pis"]])) {
        sm$parameters$pis <- sm$parameters$pis[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["means"]])) {
        sm$parameters$means <- sm$parameters$means[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["covariances"]])) {
        sm$parameters$covariances <-
          sm$parameters$covariances[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["ses"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$ses[1, 1:K] <- sm$parameters$ses[1, new_order]
        } else {
          sm$parameters$ses <- sm$parameters$ses[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["betas"]])) {
        if (length(dim(sm$parameters$betas)) == 3) {
          sm$parameters$betas <- sm$parameters$betas[new_order, , , drop = FALSE]
        } else {
          sm$parameters$betas <- sm$parameters$betas[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["beta_pooled"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$beta_pooled[1, 1:K] <- sm$parameters$beta_pooled[1, new_order]
        } else {
          # Guard: when beta_pooled is a degenerate 0x0 placeholder (produced by
          # a constant-outcome distal_pooled model), there is nothing to reorder.
          if (nrow(sm$parameters$beta_pooled) > 0 &&
              ncol(sm$parameters$beta_pooled) >= K) {
            sm$parameters$beta_pooled[, 1:K] <-
              sm$parameters$beta_pooled[, new_order, drop = FALSE]
          }
        }
      }
      return(sm)
    }

    if (inherits(model_state$sm, "nested")) {
      for (name in names(model_state$sm$models))
        model_state$sm$models[[name]] <- sort_sm_params(model_state$sm$models[[name]])
    } else {
      model_state$sm <- sort_sm_params(model_state$sm)
    }
  }

  if (!is.null(model_state$log_resp))
    model_state$log_resp <- model_state$log_resp[, new_order, drop = FALSE]

  return(model_state)
}

#' Print Measurement Model Parameters
#'
#' @description
#' Prints a formatted table of the fitted measurement model parameters:
#' item-response probabilities for categorical models, or means for Gaussian
#' models. Results are broken down by latent class. Handles both flat and
#' nested (mixed) measurement models.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#'
#' @return Invisibly returns \code{NULL}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' measurement_summary(fit)
#'
#' @export
measurement_summary <- function(object) {
  K <- object$n_components
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")

  print_item_matrix <- function(mat, title, sub_model = NULL) {
    cat(sprintf("\n%s\n", title))
    item_names <- colnames(mat)

    if (is.null(item_names)) {
      if (!is.null(sub_model) && !is.null(sub_model$max_val)) {
        M <- sub_model$max_val
        n_items <- ncol(mat) / M
        item_names <- paste0("Poly_Item_", rep(1:n_items, each = M),
                             " (Cat ", 1:M, ")")
      } else {
        item_names <- paste0("Item_", 1:ncol(mat))
      }
    }

    cat(sprintf("%-20s", "Indicator"))
    for (k in 1:K) cat(sprintf(" | Class %d", k))
    cat("\n")
    cat(paste0(rep("-", 20 + K * 10), collapse = ""), "\n")

    for (j in 1:ncol(mat)) {
      cat(sprintf("%-20s", item_names[j]))
      for (k in 1:K) cat(sprintf(" | %7.3f", mat[k, j]))
      cat("\n")
    }
  }

  mm <- object$mm
  if (inherits(mm, "nested")) {
    for (name in names(mm$models)) {
      sub_mm <- mm$models[[name]]
      if (!is.null(sub_mm$parameters$pis))
        print_item_matrix(sub_mm$parameters$pis,
                          paste("Categorical Probabilities:", toupper(name)), sub_mm)
      if (!is.null(sub_mm$parameters$means))
        print_item_matrix(sub_mm$parameters$means,
                          paste("Continuous Means:", toupper(name)), sub_mm)
    }
  } else {
    if (!is.null(mm$parameters$pis))
      print_item_matrix(mm$parameters$pis, "CATEGORICAL PROBABILITIES", mm)
    if (!is.null(mm$parameters$means))
      print_item_matrix(mm$parameters$means, "CONTINUOUS MEANS", mm)
  }
  cat("=========================================================\n")
}

#' Print Classification Diagnostics
#'
#' @description
#' Computes and prints the Average Posterior Probability (AvePP) matrix.
#' Each row corresponds to observations modally assigned to a given class;
#' each column shows the mean posterior probability for that class. High
#' values on the diagonal (and low values off it) indicate well-separated,
#' clearly-assigned classes.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#'
#' @return A K x K numeric matrix of average posterior probabilities,
#'   returned invisibly. The matrix is also printed to the console.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' classification_diagnostics(fit)
#'
#' @export
classification_diagnostics <- function(object) {
  resp       <- exp(object$log_resp)
  pred_class <- max.col(resp)
  K          <- object$n_components

  ave_pp <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    idx <- which(pred_class == k)
    if (length(idx) > 0)
      ave_pp[k, ] <- colMeans(resp[idx, , drop = FALSE])
  }
  rownames(ave_pp) <- paste("Assigned Class", 1:K)
  colnames(ave_pp) <- paste("Prob C", 1:K)

  cat("=========================================================\n")
  cat("          AVERAGE POSTERIOR PROBABILITIES (AvePP)        \n")
  cat("=========================================================\n")
  cat("Rows: Modal Assignment | Columns: Mean Probability\n\n")
  print(round(ave_pp, 3))
  cat("=========================================================\n")
  invisible(ave_pp)
}

#' Summarise a Fitted Mixture Model
#'
#' @description
#' Prints a detailed summary of the structural model parameters. Depending on
#' which structural model was fitted, this includes covariate regression
#' coefficients (as odds ratios with 95\% confidence intervals and p-values),
#' distal outcome means, or class-specific regression effects. If no structural
#' model is present, a notice is printed directing the user to
#' \code{\link{measurement_summary}}.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ref_class Integer. The reference latent class for pairwise contrasts.
#'   Defaults to the first class (\code{1}).
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly returns \code{NULL}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' summary(fit)
#'
#' @export
summary.mixture_model <- function(object, ref_class = NULL, ...) {
  K <- object$n_components
  if (is.null(ref_class)) ref_class <- 1

  # Input validation: ref_class must be a valid class index.
  # Without this guard the function starts printing output, then crashes
  # mid-way with a cryptic "subscript out of bounds" error.
  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  if (is.null(object$sm)) {
    cat("Notice: No structural model found. Use measurement_summary() for item parameters.\n")
    return(invisible())
  }

  cat("=========================================================\n")
  cat("             STRUCTURAL MODEL SUMMARY                    \n")
  cat("=========================================================\n")

  # A. Covariate model
  sm_sub <- NULL
  if (inherits(object$sm, "covariate")) sm_sub <- object$sm
  if (inherits(object$sm, "nested") && "predictor" %in% names(object$sm$models))
    sm_sub <- object$sm$models$predictor

  if (!is.null(sm_sub) && !is.null(sm_sub$parameters$hessian)) {
    cat("\nCATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)\n")
    cat(sprintf("Reference Class: %d\n", ref_class))
    cat("---------------------------------------------------------\n")
    cat("                     OR       [95% CI]        P-Value\n")

    betas     <- sm_sub$parameters$beta
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else paste0("V", 1:D)
    Sigma     <- pinv(-sm_sub$parameters$hessian)

    for (c in setdiff(1:K, ref_class)) {
      cat(sprintf("\nClass %d ON\n", c))
      for (v in 1:D) {
        est     <- betas[c, v] - betas[ref_class, v]
        idx_c   <- (c - 1) * D + v
        idx_ref <- (ref_class - 1) * D + v
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
          2 * Sigma[idx_c, idx_ref]
        se    <- sqrt(max(0, var_diff))
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("  %-15s %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                    var_names[v], exp(est),
                    exp(est - 1.96 * se), exp(est + 1.96 * se), p_val))
      }
    }
  }

  # B. Distal pooled
  pooled_sub <- NULL
  if (inherits(object$sm, "distal_pooled")) pooled_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_pooled"))
    pooled_sub <- object$sm$models$distal

  if (!is.null(pooled_sub) && !is.null(pooled_sub$parameters$beta_pooled)) {
    cat("\nPOOLED DISTAL REGRESSION (MAIN EFFECTS)\n")
    cat(sprintf("Reference Class: %d\n", ref_class))
    cat("---------------------------------------------------------\n")

    beta_mat   <- pooled_sub$parameters$beta_pooled
    M_minus_1  <- nrow(beta_mat)
    K_distal   <- K
    D_cov      <- ncol(beta_mat) - K_distal
    Sigma      <- pinv(-pooled_sub$parameters$hessian)
    var_names  <- if (D_cov > 0) paste0("Z", 1:D_cov) else character(0)

    for (m in 1:M_minus_1) {
      cat(sprintf("\nOutcome Category %d (vs Category 1) ON\n", m + 1))
      cat("                     OR       [95% CI]        P-Value\n")
      cat("\n  Latent Class (Main Effect):\n")
      for (c in setdiff(1:K_distal, ref_class)) {
        est     <- beta_mat[m, c] - beta_mat[m, ref_class]
        idx_c   <- (m - 1) * (K_distal + D_cov) + c
        idx_ref <- (m - 1) * (K_distal + D_cov) + ref_class
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
          2 * Sigma[idx_c, idx_ref]
        se    <- sqrt(max(0, var_diff))
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("    Class %d        %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                    c, exp(est), exp(est - 1.96*se), exp(est + 1.96*se), p_val))
      }
      if (D_cov > 0) {
        cat("\n  Covariates (Main Effect):\n")
        for (v in 1:D_cov) {
          est  <- beta_mat[m, K_distal + v]
          idx  <- (m - 1) * (K_distal + D_cov) + K_distal + v
          se   <- sqrt(max(0, Sigma[idx, idx]))
          z_val <- est / se
          p_val <- 2 * (1 - pnorm(abs(z_val)))
          cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                      var_names[v], exp(est),
                      exp(est - 1.96*se), exp(est + 1.96*se), p_val))
        }
      }
    }
  }

  # C. Distal regression (moderated)
  distal_sub <- NULL
  if (inherits(object$sm, "distal_regression")) distal_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_regression"))
    distal_sub <- object$sm$models$distal

  if (!is.null(distal_sub) && !is.null(distal_sub$parameters$betas)) {
    cat("\nSIMULTANEOUS DISTAL REGRESSION (MODERATED BY CLASS)\n")
    cat("---------------------------------------------------------\n")

    distal_betas <- distal_sub$parameters$betas
    K_distal     <- dim(distal_betas)[1]
    M_minus_1    <- dim(distal_betas)[2]
    D_distal     <- dim(distal_betas)[3]

    # When the outcome is constant (M=1), init_params stores a
    # degenerate array with dim c(K, 0, 0), so M_minus_1 = D_distal = 0.
    # In R, 1:0 evaluates to c(1, 0) — NOT an empty sequence — so the inner
    # loops fire on a zero-dimension array and crash with "subscript out of
    # bounds".  seq_len(0) correctly returns integer(0), skipping the loops.
    if (M_minus_1 == 0) {
      cat("  (Constant outcome - no parameters to display)\n")
    } else {
      var_names <- c("Intercept", paste0("Z", seq_len(D_distal - 1)))

      for (k in seq_len(K_distal)) {
        cat(sprintf("\nClass %d:\n", k))
        cat("                     OR       [95% CI]        P-Value\n")
        if (!is.null(distal_sub$parameters$hessians) &&
            length(distal_sub$parameters$hessians) >= k) {
          Sigma <- pinv(-distal_sub$parameters$hessians[[k]])
        } else {
          Sigma <- matrix(0, M_minus_1 * D_distal, M_minus_1 * D_distal)
        }
        for (m in seq_len(M_minus_1)) {
          cat(sprintf("  Outcome Category %d (vs Category 1) ON\n", m + 1))
          for (v in seq_len(D_distal)) {
            est  <- distal_betas[k, m, v]
            idx  <- (m - 1) * D_distal + v
            se   <- sqrt(max(0, Sigma[idx, idx]))
            if (se > 0) {
              z_val <- est / se
              p_val <- 2 * (1 - pnorm(abs(z_val)))
              cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                          var_names[v], exp(est),
                          exp(est - 1.96*se), exp(est + 1.96*se), p_val))
            } else {
              cat(sprintf("    %-13s %7.3f  [   N/A,    N/A]       N/A\n",
                          var_names[v], exp(est)))
            }
          }
        }
      }
    }
  }

  # D. Continuous distal (means)
  cont_sub <- NULL
  if (inherits(object$sm, "distal_continuous")) cont_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous"))
    cont_sub <- object$sm$models$distal

  if (!is.null(cont_sub) && !is.null(cont_sub$parameters$means)) {
    cat("\nCONTINUOUS DISTAL OUTCOME (MEANS)\n")
    cat("---------------------------------------------------------\n")
    cat("                 Mean       [95% CI]        SE\n")
    for (k in 1:K) {
      mu <- cont_sub$parameters$means[k]
      se <- cont_sub$parameters$ses[k]
      cat(sprintf("  Class %d      %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                  k, mu, mu - 1.96 * se, mu + 1.96 * se, se))
    }
  }

  # E. Continuous distal regression
  cont_reg_sub <- NULL
  if (inherits(object$sm, "distal_continuous_regression")) cont_reg_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_regression"))
    cont_reg_sub <- object$sm$models$distal

  if (!is.null(cont_reg_sub) && !is.null(cont_reg_sub$parameters$betas)) {
    cat("\nCONTINUOUS DISTAL REGRESSION (Y ~ Z * Class)\n")
    cat("---------------------------------------------------------\n")
    betas     <- cont_reg_sub$parameters$betas
    ses       <- cont_reg_sub$parameters$ses
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else
      c("Intercept", paste0("Z", 1:(D - 1)))
    for (k in 1:K) {
      cat(sprintf("\nClass %d:\n", k))
      cat("                 Estimate   [95% CI]        P-Value\n")
      for (v in 1:D) {
        est   <- betas[k, v]
        se    <- ses[k, v]
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("  %-13s %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                    var_names[v], est, est - 1.96*se, est + 1.96*se, p_val))
      }
    }
  }

  # F. Continuous distal pooled
  cont_pool_sub <- NULL
  if (inherits(object$sm, "distal_continuous_pooled")) cont_pool_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_pooled"))
    cont_pool_sub <- object$sm$models$distal

  if (!is.null(cont_pool_sub) && !is.null(cont_pool_sub$parameters$beta_pooled)) {
    cat("\nCONTINUOUS DISTAL POOLED REGRESSION (Main Effects)\n")
    cat("---------------------------------------------------------\n")
    theta <- as.vector(cont_pool_sub$parameters$beta_pooled)
    ses   <- as.vector(cont_pool_sub$parameters$ses)
    K_distal <- K
    D_cov <- length(theta) - K_distal
    var_names <- if (D_cov > 0) paste0("Z", 1:D_cov) else character(0)

    cat("\n  Latent Class (Intercepts):\n")
    cat("                 Estimate   [95% CI]        P-Value\n")
    for (k in 1:K_distal) {
      est <- theta[k]
      se  <- ses[k]
      z_val <- est / se
      p_val <- 2 * (1 - pnorm(abs(z_val)))
      cat(sprintf("    Class %d      %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                  k, est, est - 1.96*se, est + 1.96*se, p_val))
    }

    if (D_cov > 0) {
      cat("\n  Covariates (Pooled Slopes):\n")
      for (v in 1:D_cov) {
        est <- theta[K_distal + v]
        se  <- ses[K_distal + v]
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("    %-11s %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                    var_names[v], est, est - 1.96*se, est + 1.96*se, p_val))
      }
    }
  }

  cat("=========================================================\n")
}

#' Fit a Latent Mixture Model (LCA / LPA)
#'
#' @description
#' The core estimation function. Fits a latent class analysis (LCA) or latent
#' profile analysis (LPA) model using the EM algorithm. Optionally fits a
#' structural model (covariates or distal outcomes) using 1-, 2-, or 3-step
#' estimation with optional bias correction.
#'
#' @param X A numeric matrix or data frame of indicator variables for the
#'   measurement model. Rows are observations; columns are items or variables.
#' @param Y Optional numeric matrix or data frame of outcome or covariate
#'   variables for the structural model. Must be provided when
#'   \code{structural} is not \code{NULL}.
#' @param n_components Positive integer. Number of latent classes (or profiles)
#'   to estimate. Default is \code{2}.
#' @param measurement Character string or named list specifying the measurement
#'   model type. Accepted strings: \code{"binary"} / \code{"bernoulli"},
#'   \code{"binary_nan"} / \code{"bernoulli_nan"} (handles \code{NA}s),
#'   \code{"categorical"} / \code{"multinoulli"},
#'   \code{"categorical_nan"} / \code{"multinoulli_nan"},
#'   \code{"continuous"} / \code{"gaussian_diag"},
#'   \code{"continuous_nan"} / \code{"gaussian_diag_nan"},
#'   \code{"gaussian"} / \code{"gaussian_unit"}.
#'   Pass a named list to specify a mixed (nested) measurement model with
#'   different variable types. Default is \code{"binary"}.
#' @param structural Character string specifying the structural model type.
#'   One of \code{"covariate"}, \code{"distal_regression"},
#'   \code{"distal_pooled"}, \code{"distal_continuous"},
#'   \code{"distal_continuous_regression"}. Requires \code{Y}. Default is
#'   \code{NULL} (measurement model only).
#' @param n_steps Integer. Estimation approach: \code{1} for simultaneous
#'   1-step, \code{2} for 2-step, or \code{3} for bias-corrected 3-step.
#'   Default is \code{1}.
#' @param correction Character. Bias correction for 3-step estimation.
#'   One of \code{"none"}, \code{"BCH"}, or \code{"ML"}. Ignored when
#'   \code{n_steps} is not \code{3}. Default is \code{"none"}.
#' @param n_init Positive integer. Number of random restarts. The solution
#'   with the highest log-likelihood is retained. Default is \code{1}.
#' @param max_iter Positive integer. Maximum EM iterations per restart.
#'   Default is \code{1000}.
#' @param random_state Optional integer seed for reproducibility. Default is
#'   \code{NULL}.
#' @param order_by_size Logical. If \code{TRUE} (default), classes are sorted
#'   from largest to smallest after fitting.
#' @param weights Optional numeric vector of length \code{nrow(X)} for survey
#'   or case weights. Default is \code{NULL} (equal weights of 1).
#' @param ... Additional arguments passed to the measurement or structural
#'   model constructors (e.g., \code{max_val} for multinoulli models).
#'
#' @return An object of class `mixture_model`, a list with:
#'   * `n_components` Number of latent classes.
#'   * `weights` Numeric vector of estimated class proportions.
#'   * `mm` Fitted measurement model state object.
#'   * `sm` Fitted structural model state object, or `NULL`.
#'   * `metrics` Named list: `ll` (log-likelihood), `aic`, `bic`, `sabic`,
#'     `n_params`, and `entropy` (relative entropy, 0-1 scale).
#'   * `log_resp` Matrix of log posterior class probabilities (n x K).
#'     Use `exp(fit$log_resp)` to obtain posterior probabilities.
#'   * `converged` Logical. Whether the EM algorithm converged.
#'   * `n_iter` Integer. Number of EM iterations run.
#'   * `step1_metrics` Named list of Step-1 fit indices (only when
#'     `n_steps = 3`).
#'
#' @examples
#' # Binary LCA with 3 classes and 5 random restarts
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 3, measurement = "binary", n_init = 5)
#' print(fit)
#' summary(fit)
#' measurement_summary(fit)
#'
#' # Continuous LPA (2 classes)
#' X_cont <- matrix(rnorm(300), nrow = 100)
#' fit_lpa <- fit_mixture(X_cont, n_components = 2, measurement = "continuous")
#'
#' # 3-step LCA with a covariate and ML correction
#' Z <- matrix(rnorm(100), nrow = 100)
#' fit_cov <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                        structural = "covariate",
#'                        n_steps = 3, correction = "ML", n_init = 5)
#' summary(fit_cov)
#'
#' @export
#' @importFrom stats complete.cases cov dnorm optim pchisq pnorm qnorm rbinom rnorm runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
fit_mixture <- function(X, Y = NULL, n_components = 2,
                        measurement = "binary", structural = NULL,
                        n_steps = 1, correction = "none", n_init = 1,
                        max_iter = 1000, random_state = NULL,
                        order_by_size = TRUE, weights = NULL, ...) {

  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.null(Y) && is.data.frame(Y)) Y <- as.matrix(Y)

  n_samples <- nrow(X)

  # --- Input Validation ---

  # n_components must be a positive integer (Bug 2 fix: k=0 silently produced
  # garbage; k<0 crashed deep in initialization with a cryptic message)
  if (!is.numeric(n_components) || length(n_components) != 1 ||
      !is.finite(n_components) || n_components < 1L)
    stop(sprintf(
      "n_components must be a positive integer (>= 1). Got: %s",
      n_components
    ))

  # n_steps must be 1, 2, or 3
  if (!n_steps %in% c(1L, 2L, 3L))
    stop(sprintf("n_steps must be 1, 2, or 3. Got: %d", n_steps))

  # correction must be one of the three supported values (Gap 10 fix)
  valid_corrections <- c("none", "ML", "BCH")
  if (!correction %in% valid_corrections)
    stop(sprintf(
      "correction '%s' not recognized. Choose from: %s",
      correction, paste(valid_corrections, collapse = ", ")
    ))

  # Warn when X contains NAs but a non-NaN measurement family is used (Bug 7 fix)
  if (anyNA(X) && is.character(measurement)) {
    nan_variants <- c("bernoulli_nan", "binary_nan",
                      "multinoulli_nan", "categorical_nan",
                      "gaussian_diag_nan", "continuous_nan")
    if (!measurement %in% nan_variants)
      warning(sprintf(
        paste0("X contains NA values but measurement = '%s' does not handle missing data. ",
               "Did you mean '%s_nan'? Proceeding, but results may contain NAs."),
        measurement, measurement
      ))
  }

  # Validate binary data when a Bernoulli family is requested (Bug 8 fix)
  if (is.character(measurement) &&
      measurement %in% c("binary", "bernoulli", "binary_nan", "bernoulli_nan")) {
    valid_vals <- X[!is.na(X)]
    if (length(valid_vals) > 0 && !all(valid_vals %in% c(0, 1)))
      stop(sprintf(
        paste0("measurement = '%s' requires X values in {0, 1}. ",
               "Found values outside this set: %s"),
        measurement,
        paste(sort(unique(valid_vals[!valid_vals %in% c(0, 1)]))[1:min(5, sum(!valid_vals %in% c(0,1)))],
              collapse = ", ")
      ))
  }

  if (is.null(weights)) {
    weights <- rep(1, n_samples)
  } else {
    if (length(weights) != n_samples)
      stop("Length of weights must match rows of X.")
  }

  # Structural model requires Y. Without this guard the SM is built but never
  # fitted (m_step_core gates on !is.null(Y)), so parameters$beta stays NULL
  # and every downstream function (coef, confint, wald tests) crashes with a
  # cryptic error rather than pointing here.
  if (!is.null(structural) && is.null(Y))
    stop(paste(
      "A structural model was specified but Y is NULL.",
      "Provide a Y matrix containing the outcome/covariate data,",
      "or set structural = NULL for a measurement-only model."
    ))

  model_state <- list(
    n_components          = n_components,
    weights               = rep(1 / n_components, n_components),
    mm                    = build_emission(measurement, n_components = n_components, ...),
    sm                    = if (!is.null(structural))
      build_emission(structural, n_components = n_components, ...)
    else NULL,
    n_steps               = n_steps,
    correction            = correction,
    sample_weights        = weights,
    # Store the original descriptor so bootstrap.R can re-fit replicates
    # using the same measurement specification (Bug 2 fix).
    measurement_descriptor = measurement
  )
  class(model_state) <- "mixture_model"

  if (n_steps == 1) {
    model_state <- fit_em(model_state, X, Y, n_init, max_iter, random_state)

  } else if (n_steps == 2) {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state)
    if (!is.null(Y) && !is.null(model_state$sm)) {
      resp <- exp(model_state$log_resp)
      model_state$sm <- init_params(model_state$sm, Y, resp)
      model_state$sm <- m_step(model_state$sm, Y, resp)
    }

  } else if (n_steps == 3) {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state)

    # Step 1 metrics (measurement model only)
    n_params_s1  <- n_parameters(model_state$mm) + (model_state$n_components - 1)
    ll_s1        <- sum(model_state$sample_weights * model_state$lower_bound)
    resp_s1      <- exp(model_state$log_resp)
    abs_ent_s1   <- sum(model_state$sample_weights *
                          (-resp_s1 * log(resp_s1 + 1e-15)))
    # Use relative_entropy() helper to handle K=1 correctly (Gap 9 fix)
    rel_ent_s1   <- relative_entropy(abs_ent_s1,
                                     sum(model_state$sample_weights),
                                     model_state$n_components)

    model_state$step1_metrics <- list(
      ll       = ll_s1,
      n_params = n_params_s1,
      aic      = -2 * ll_s1 + 2 * n_params_s1,
      bic      = -2 * ll_s1 + log(n_samples) * n_params_s1,
      sabic    = -2 * ll_s1 + log((n_samples + 2) / 24) * n_params_s1,
      entropy  = rel_ent_s1
    )

    if (!is.null(Y) && !is.null(model_state$sm)) {
      if (correction == "ML") {
        model_state <- fit_ml(model_state, X, Y, max_iter = max_iter)
      } else if (correction == "BCH") {
        model_state <- fit_bch(model_state, X, Y)
      } else {
        # correction = "none": plain 2-step update on the structural model (Bug 1 fix).
        # The measurement model is already frozen; we just fit the SM on the
        # posterior responsibilities from step 1 without any bias correction.
        resp <- exp(model_state$log_resp)
        model_state$sm <- init_params(model_state$sm, Y, resp)
        model_state$sm <- m_step(model_state$sm, Y, resp)
      }
    }
  }

  if (order_by_size) model_state <- sort_model_classes(model_state)

  # Attach column names.
  # Guard: only assign colnames to pis when dimensions match (Bernoulli: J cols;
  # Multinoulli: J*M cols — colnames(X) has length J so would misfire there).
  if (!is.null(colnames(X)) && !is.null(model_state$mm$parameters$pis) &&
      ncol(model_state$mm$parameters$pis) == ncol(X))
    colnames(model_state$mm$parameters$pis) <- colnames(X)

  # Attach covariate names to beta only when the SM is initialized (Bug 1 fix:
  # guard against NULL beta after correction="none" before the else branch above
  # ran, or any other path where the SM was never fit).
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "covariate") &&
      !is.null(model_state$sm$parameters$beta)) {
    intercept_flag <- isTRUE(model_state$sm$intercept)
    expected_D     <- ncol(model_state$sm$parameters$beta)

    y_col_names <- if (!is.null(colnames(Y))) colnames(Y)
    else paste0("V", seq_len(ncol(Y)))
    cov_names   <- if (intercept_flag) c("Intercept", y_col_names)
    else y_col_names
    if (!is.null(cov_names) && !is.null(expected_D) &&
        length(cov_names) == expected_D)
      colnames(model_state$sm$parameters$beta) <- cov_names
  }

  # Final metrics
  n_params <- n_parameters(model_state$mm) + (model_state$n_components - 1)
  if (!is.null(model_state$sm)) n_params <- n_params + n_parameters(model_state$sm)
  ll       <- sum(model_state$sample_weights * model_state$lower_bound)
  resp     <- exp(model_state$log_resp)
  abs_ent  <- sum(model_state$sample_weights * (-resp * log(resp + 1e-15)))
  # Use relative_entropy() helper to handle K=1 (returns 1) cleanly (Gap 9 fix)
  ent      <- relative_entropy(abs_ent,
                               sum(model_state$sample_weights),
                               model_state$n_components)

  model_state$metrics <- list(
    ll       = ll,
    n_params = n_params,
    aic      = -2 * ll + 2 * n_params,
    bic      = -2 * ll + log(n_samples) * n_params,
    sabic    = -2 * ll + log((n_samples + 2) / 24) * n_params,
    entropy  = ent
  )

  return(model_state)
}

#' Print a Brief Summary of a Fitted Mixture Model
#'
#' @description
#' Prints a compact overview of the fitted model including: number of classes,
#' estimation method, convergence status, log-likelihood, relative entropy,
#' and estimated class proportions. For full parameter tables, use
#' \code{\link{summary.mixture_model}} (structural parameters) or
#' \code{\link{measurement_summary}} (item parameters).
#'
#' @param x A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly returns \code{x}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' print(fit)
#' # or equivalently:
#' fit
#'
#' @export
print.mixture_model <- function(x, ...) {
  cat("=========================================================\n")
  cat("                  LATENT MIXTURE MODEL                   \n")
  cat("=========================================================\n")
  cat(sprintf("Classes Estimated  : %d\n", x$n_components))
  cat(sprintf("Estimation Method  : %d-step\n", x$n_steps))
  if (x$n_steps == 3) cat(sprintf("Correction Method  : %s\n", x$correction))
  cat(sprintf("Converged          : %s (in %d iterations)\n", x$converged, x$n_iter))
  cat("---------------------------------------------------------\n")
  if (!is.null(x$step1_metrics)) {
    cat(sprintf("  Log-Likelihood (Step 1) : %.2f\n", x$step1_metrics$ll))
    cat(sprintf("  Rel. Entropy   (Step 1) : %.4f\n", x$step1_metrics$entropy))
  } else {
    cat(sprintf("  Log-Likelihood : %.2f\n", x$metrics$ll))
    cat(sprintf("  Rel. Entropy   : %.4f\n", x$metrics$entropy))
  }
  cat("---------------------------------------------------------\n")
  cat("Class Weights (Sizes):\n")
  for (i in seq_along(x$weights))
    cat(sprintf("  Class %d: %.2f%%\n", i, x$weights[i] * 100))
  cat("=========================================================\n")
  cat("Type summary(model) for structural parameters or measurement_summary(model) for item parameters.\n")
}

#' Compare Mixture Models Across a Range of Class Numbers
#'
#' @description
#' Fits a sequence of measurement-only mixture models, one for each value of
#' \code{k} in \code{k_range}, and returns a table of fit indices to guide
#' class enumeration. The best model according to BIC is identified
#' automatically.
#'
#' @param X A numeric matrix or data frame of indicator variables.
#' @param k_range Integer vector of class numbers to fit. All values must be >= 1. Default is \code{1:5}.
#' @param measurement Character string specifying the measurement model type.
#'   See \code{\link{fit_mixture}} for accepted values. Default is
#'   \code{"binary"}.
#' @param n_init Positive integer. Number of random restarts per model.
#'   Default is \code{10}.
#' @param n_steps Integer. Estimation method: \code{1}, \code{2}, or \code{3}.
#'   Default is \code{1}.
#' @param ... Additional arguments passed to \code{\link{fit_mixture}}.
#'
#' @return A named list with three elements:
#'   * `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
#'     `Params`, `AIC`, `BIC`, `SABIC`, `Entropy`.
#'   * `models` Named list of fitted `mixture_model` objects, one per K
#'     (names are `"K1"`, `"K2"`, etc.).
#'   * `best_k` Integer. The value of K with the lowest BIC.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
#'                            n_init = 5)
#' result$fit_table
#' result$best_k
#'
#' @export
compare_mixtures <- function(X, k_range = 1:5, measurement = "binary",
                             n_init = 10, n_steps = 1, ...) {
  # k=0 would silently fit a degenerate model with LL=-Inf; negative
  # k values crash deep in initialisation with a cryptic error.
  if (any(k_range < 1L))
    stop(sprintf(
      "All values in k_range must be >= 1. Got invalid values: %s",
      paste(sort(unique(k_range[k_range < 1L])), collapse = ", ")
    ))
  cat(sprintf("Running Model Selection across K = %d to %d...\n\n",
              min(k_range), max(k_range)))
  results <- list()
  models  <- list()
  for (k in k_range) {
    cat(sprintf("Fitting %d-class model...\n", k))
    fit <- fit_mixture(X = X, Y = NULL, n_components = k,
                       measurement = measurement,
                       n_steps = n_steps, n_init = n_init, ...)
    models[[paste0("K", k)]] <- fit
    results[[k]] <- data.frame(
      Classes = k, LL = fit$metrics$ll, Params = fit$metrics$n_params,
      AIC = fit$metrics$aic, BIC = fit$metrics$bic,
      SABIC = fit$metrics$sabic, Entropy = fit$metrics$entropy
    )
  }
  fit_table   <- do.call(rbind, results)
  best_bic_k  <- fit_table$Classes[which.min(fit_table$BIC)]
  cat("\n=== Model Selection Summary ===\n")
  print(round(fit_table, 3))
  cat(sprintf("\n-> Best model according to BIC: %d classes\n", best_bic_k))
  return(list(fit_table = fit_table, models = models, best_k = best_bic_k))
}

#' Extract Covariate Odds Ratios from a Fitted Mixture Model
#'
#' @description
#' Extracts the logistic regression coefficients from a covariate structural
#' model and returns them as a matrix of odds ratios, centered on a reference
#' class. Only available when the model was fitted with
#' \code{structural = "covariate"}.
#'
#' @param object A fitted \code{mixture_model} object with a covariate
#'   structural model.
#' @param ref_class Integer. The reference class for centering. All other
#'   class odds ratios are expressed relative to this class. Default is
#'   \code{1}.
#' @param covariate_names Optional character vector of predictor names to
#'   override the column names stored in the model. Default is \code{NULL}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return A K x D numeric matrix of odds ratios, where rows are latent
#'   classes and columns are predictors (including the intercept). The
#'   reference class row will always show \code{1.000}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' coef(fit)
#'
#' @export
coef.mixture_model <- function(object, ref_class = 1, covariate_names = NULL, ...) {
  if (is.null(object$sm) || !inherits(object$sm, "covariate"))
    stop("No covariate model.")
  K     <- object$n_components
  betas <- object$sm$parameters$beta
  if (!is.null(covariate_names))
    colnames(betas) <- c("Intercept", covariate_names)
  betas_ref   <- sweep(betas, 2, betas[ref_class, ], "-")
  odds_ratios <- exp(betas_ref)
  rownames(odds_ratios) <- paste("Class", 1:K)
  rownames(odds_ratios)[ref_class] <- paste("Class", ref_class, "(Ref)")
  return(odds_ratios)
}
