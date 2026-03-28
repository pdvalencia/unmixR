# ==============================================================================
# unmixR Implementation - Utilities and Math Helpers
# ==============================================================================

# Log-Sum-Exp trick for numerical stability
logsumexp <- function(x, MARGIN = 1) {
  if (is.null(dim(x))) {
    max_x <- max(x)
    if (max_x == -Inf) return(-Inf)
    return(max_x + log(sum(exp(x - max_x))))
  } else {
    max_x <- apply(x, MARGIN, max)
    max_x[max_x == -Inf] <- 0
    if (MARGIN == 1) {
      res <- max_x + log(rowSums(exp(x - max_x)))
    } else {
      res <- max_x + log(colSums(exp(sweep(x, 2, max_x, "-"))))
    }
    return(res)
  }
}

# Input validation: Check if positive
check_positive <- function(...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!is.numeric(val) || any(val <= 0, na.rm = TRUE))
      stop(sprintf("Expected %s > 0, got %s", name, paste(val, collapse=" ")), call. = FALSE)
  }
}

# Input validation: Check if non-negative
check_nonneg <- function(...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!is.numeric(val) || any(val < 0, na.rm = TRUE))
      stop(sprintf("Expected %s >= 0, got %s", name, paste(val, collapse=" ")), call. = FALSE)
  }
}

# Input validation: Check if value is in allowed choices
check_in <- function(choices, ...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!all(val %in% choices))
      stop(sprintf("%s value '%s' not recognized. Choose from: %s",
                   name, paste(val, collapse=" "), paste(choices, collapse=", ")),
           call. = FALSE)
  }
}

# Modal (Hard) Assignment of Probabilities
modal <- function(resp, clip = FALSE) {
  max_idx <- max.col(resp, ties.method = "random")
  modal_resp <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
  modal_resp[cbind(1:nrow(resp), max_idx)] <- 1
  if (clip) modal_resp <- pmax(pmin(modal_resp, 1 - 1e-15), 1e-15)
  return(modal_resp)
}

# Multiple categorical one-hot encoding (0-indexed input)
max_one_hot <- function(mat, max_n_outcomes = NULL, total_outcomes = NULL) {
  n_samples <- nrow(mat)
  n_features <- ncol(mat)

  if (is.null(max_n_outcomes) || is.null(total_outcomes)) {
    outcomes <- apply(mat, 2, max, na.rm = TRUE) + 1
    total_outcomes <- sum(outcomes)
    max_n_outcomes <- max(outcomes)
  }

  one_hot <- matrix(0, nrow = n_samples, ncol = n_features * max_n_outcomes)

  for (c in 1:n_features) {
    integer_codes <- mat[, c]
    not_observed <- is.na(integer_codes)
    integer_codes[not_observed] <- 0
    col_indices <- integer_codes + (c - 1) * max_n_outcomes + 1
    one_hot[cbind(1:n_samples, col_indices)] <- 1
    if (any(not_observed)) {
      start_col <- (c - 1) * max_n_outcomes + 1
      end_col   <- c * max_n_outcomes
      one_hot[not_observed, start_col:end_col] <- NA
    }
  }

  return(list(one_hot = one_hot,
              max_n_outcomes = max_n_outcomes,
              total_outcomes = total_outcomes))
}

# Moore-Penrose pseudoinverse (pure base R)
pinv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X)
  d <- s$d
  d_inv <- ifelse(d > tol * d[1], 1/d, 0)
  return(s$v %*% diag(d_inv, nrow = length(d)) %*% t(s$u))
}

# Relative Entropy (normalised by log(K), scales 0-1)
relative_entropy <- function(absolute_entropy, n_samples, n_classes) {
  if (n_classes <= 1) return(1)
  rel_ent <- 1 - (absolute_entropy / (n_samples * log(n_classes)))
  return(max(0, min(1, rel_ent)))
}

# Mean imputation for missing covariates in structural models
impute_covariates <- function(Z) {
  Z <- as.matrix(Z)
  if (ncol(Z) == 0) return(Z)
  for (j in seq_len(ncol(Z))) {
    na_idx <- is.na(Z[, j])
    if (any(na_idx)) {
      col_mean <- mean(Z[, j], na.rm = TRUE)
      # All-NA column: mean() returns NaN, which would silently propagate
      # through logit calculations. Warn and impute 0 instead.
      if (is.nan(col_mean)) {
        warning(sprintf(
          paste0("impute_covariates: column %d is entirely NA. ",
                 "Imputing with 0. Check your data for completely missing covariates."),
          j
        ))
        col_mean <- 0
      }
      Z[na_idx, j] <- col_mean
    }
  }
  return(Z)
}

# ==============================================================================
# Shared helpers for categorical distal outcome models
# (used by both distal_regression.R and distal_pooled.R)
# ==============================================================================

# One-hot encode a 1-indexed integer vector Y into an (n x max_val) matrix.
# Missing values produce all-zero rows.
distal_one_hot <- function(Y, max_val) {
  n   <- length(Y)
  out <- matrix(0, nrow = n, ncol = max_val)
  valid <- !is.na(Y)
  if (any(valid)) out[cbind(which(valid), Y[valid])] <- 1
  return(out)
}

# Multinomial softmax forward pass.
# The first category is the reference (fixed logit = 0).
# @param Z        N x D design matrix
# @param beta_matrix (M-1) x D coefficient matrix
distal_forward <- function(Z, beta_matrix) {
  logits_active <- Z %*% t(beta_matrix)
  logits_full   <- cbind(0, logits_active)
  max_logits    <- apply(logits_full, 1, max)
  exp_logits    <- exp(sweep(logits_full, 1, max_logits, "-"))
  return(sweep(exp_logits, 1, rowSums(exp_logits), "/"))
}

# ==============================================================================
# prepare_covariates()
#
# Converts a user-supplied Y (data.frame, matrix, or vector) into a plain
# numeric matrix suitable for structural models, with two behaviours:
#
#   * Numeric / integer columns are passed through unchanged.
#   * Factor / character columns are dummy-coded (reference = first level):
#       - 2-level factor  → 1 binary column named after the variable
#       - k-level factor  → k-1 columns named "var.level2", "var.level3", ...
#
# Column names are always preserved / generated, so downstream summary
# functions can display real variable names instead of Z1, Z2, etc.
# ==============================================================================
prepare_covariates <- function(Y) {
  # Plain numeric matrix — nothing to do
  if (is.matrix(Y) && is.numeric(Y)) return(Y)

  # Bare numeric vector — wrap in single-column matrix
  if (is.numeric(Y) && is.null(dim(Y)))
    return(matrix(Y, ncol = 1L,
                  dimnames = list(NULL, deparse(substitute(Y)))))

  df <- as.data.frame(Y)
  out_cols <- vector("list", ncol(df))
  idx <- 0L

  for (nm in names(df)) {
    col <- df[[nm]]

    # ── numeric / integer: pass through ──────────────────────────────────────
    if (is.numeric(col) || is.integer(col)) {
      idx <- idx + 1L
      m   <- matrix(as.numeric(col), ncol = 1L,
                    dimnames = list(NULL, nm))
      out_cols[[idx]] <- m
      next
    }

    # ── factor / character: dummy-code ───────────────────────────────────────
    if (!is.factor(col)) col <- factor(col)
    lvls <- levels(col)

    if (length(lvls) < 2L) {
      warning(sprintf(
        "prepare_covariates: variable '%s' has fewer than 2 levels and will be dropped.",
        nm))
      next
    }

    other <- lvls[-1L]   # reference = first level

    if (length(other) == 1L) {
      # Binary: single column named "var.Level" so the active level is clear
      idx <- idx + 1L
      out_cols[[idx]] <- matrix(as.integer(col == other),
                                ncol = 1L,
                                dimnames = list(NULL, paste0(nm, ".", other)))
    } else {
      # Multicategorical: k-1 columns named "var.levelX"
      idx <- idx + 1L
      m   <- matrix(0L, nrow = length(col), ncol = length(other),
                    dimnames = list(NULL, paste0(nm, ".", other)))
      for (j in seq_along(other))
        m[, j] <- as.integer(col == other[j])
      out_cols[[idx]] <- m
    }
  }

  do.call(cbind, out_cols[seq_len(idx)])
}
