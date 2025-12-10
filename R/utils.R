
################################################################################
create_stratified_folds <- function(y, k) {
  folds <- vector("list", k)

  # For each class/level
  for (class in unique(y)) {
    class_indices <- which(y == class)
    n_class <- length(class_indices)

    # Shuffle and assign to folds
    shuffled <- sample(class_indices)
    fold_assignment <- rep(1:k, length.out = n_class)

    # Add to folds
    for (fold_id in 1:k) {
      folds[[fold_id]] <- c(folds[[fold_id]], shuffled[fold_assignment == fold_id])
    }
  }

  # Shuffle within each fold (optional, for randomness)
  folds <- lapply(folds, sample)

  return(folds)
}


# y <- factor(c(rep("A", 100), rep("B", 200), rep("C", 50)))
# folds <- create_stratified_folds(y, k = 5)
# lapply(folds, function(idx) table(y[idx]))


################################################################################

extract_rapid_metrics <- function(res,
                                  model,
                                  sensitive,
                                  eval_method,
                                  type = c("categorical", "numeric")) {

  type <- match.arg(type)

  # -----------------------------
  # NUMERIC
  # -----------------------------
  if (type == "numeric") {

    return(data.frame(
      model              = model,
      sensitive_variable = sensitive,
      eval_method        = "numeric",

      confidence_rate = res$confidence_rate,
      n_at_risk        = res$n_at_risk,
      percentage       = res$percentage,

      normalized_gain  = NA_real_,
      normalized_ce    = NA_real_,
      baseline         = NA_real_,
      relative_score   = NA_real_,

      stringsAsFactors = FALSE
    ))
  }

  # -----------------------------
  # CATEGORICAL
  # -----------------------------
  method <- res$method

  # --- RCS_conditional ---
  if (method == "RCS_conditional") {

    return(data.frame(
      model              = model,
      sensitive_variable = sensitive,
      eval_method        = eval_method,

      confidence_rate = res$confidence_rate,
      n_at_risk        = res$n_at_risk,
      percentage       = res$percentage,

      relative_score   = res$relative_score,
      baseline         = res$baseline,

      normalized_gain  = NA_real_,
      normalized_ce    = NA_real_,

      stringsAsFactors = FALSE
    ))
  }

  # --- RCS_marginal ---
  if (method == "RCS_marginal") {

    return(data.frame(
      model              = model,
      sensitive_variable = sensitive,
      eval_method        = eval_method,

      confidence_rate = res$confidence_rate,
      n_at_risk        = res$n_at_risk,
      percentage       = res$percentage,

      normalized_gain  = res$normalized_gain,
      baseline         = res$baseline,

      relative_score   = NA_real_,
      normalized_ce    = NA_real_,

      stringsAsFactors = FALSE
    ))
  }

  # --- NCE ---
  if (method == "NCE") {

    return(data.frame(
      model              = model,
      sensitive_variable = sensitive,
      eval_method        = eval_method,

      confidence_rate = res$risk_rate,
      n_at_risk        = NA_integer_,
      percentage       = 100 * res$risk_rate,

      normalized_ce    = res$normalized_ce,

      relative_score   = NA_real_,
      normalized_gain  = NA_real_,
      baseline         = NA_real_,

      stringsAsFactors = FALSE
    ))
  }

  stop("Unknown categorical evaluation method.")
}

