# ==============================================================================
# S3 Emission Factory (List Dispatcher)
# ==============================================================================

# Build an emission model state based on string or list descriptor
build_emission <- function(descriptor, n_components = 2, ...) {

  # 1. Handle Nested Models (List of sub-models)
  if (is.list(descriptor)) {
    state <- list(
      n_components = n_components,
      models = list(),
      columns_per_model = numeric()
    )
    class(state) <- c("nested", "emission")

    for (name in names(descriptor)) {
      item <- descriptor[[name]]
      model_type <- item$model
      n_cols <- item$n_columns

      args <- item
      args$model <- NULL
      args$n_columns <- NULL
      args$n_components <- n_components

      # Build the sub-model state recursively
      state$models[[name]] <- do.call(build_emission, c(list(descriptor = model_type), args))
      state$columns_per_model[name] <- n_cols
    }
    return(state)
  }

  # 2. Handle Simple Models (Character string)
  if (is.character(descriptor)) {
    args_list <- list(n_components = n_components, ...)

    if (descriptor %in% c("bernoulli", "binary")) {
      return(do.call(categorical_model, c(list(type = "bernoulli"), args_list)))
    } else if (descriptor %in% c("bernoulli_nan", "binary_nan")) {
      return(do.call(categorical_model, c(list(type = "bernoulli_nan"), args_list)))
    } else if (descriptor %in% c("multinoulli", "categorical")) {
      return(do.call(categorical_model, c(list(type = "multinoulli"), args_list)))
    } else if (descriptor %in% c("multinoulli_nan", "categorical_nan")) {
      return(do.call(categorical_model, c(list(type = "multinoulli_nan"), args_list)))
    } else if (descriptor %in% c("gaussian_unit", "gaussian")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_unit"), args_list)))
    } else if (descriptor %in% c("gaussian_diag", "continuous")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_diag"), args_list)))
    } else if (descriptor %in% c("gaussian_diag_nan", "continuous_nan")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_diag_nan"), args_list)))
    } else if (descriptor == "covariate") {
      return(do.call(covariate_model, args_list))
    } else if (descriptor == "distal_regression") {
      return(do.call(distal_regression_model, args_list))
    } else if (descriptor == "distal_pooled") {
      return(do.call(distal_pooled_model, args_list))
    } else if (descriptor == "distal_continuous") {
      return(do.call(distal_continuous_model, args_list))
    } else if (descriptor == "distal_continuous_regression") {
      return(do.call(distal_continuous_regression_model, args_list))
    } else {
      stop(sprintf("Emission descriptor '%s' not recognized.", descriptor))
    }
  }
}
