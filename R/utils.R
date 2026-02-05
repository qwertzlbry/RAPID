#' @importFrom graphics grid points text
#' @importFrom stats binomial complete.cases glm median model.matrix qnorm sd
NULL
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
################################################################################
#' Print method for RAPID results
#'
#' @param x A rapid_result object
#' @param type Type of output: "summary" (default) or "high_risk"
#' @param ... Additional arguments (unused)
#' @export
#'
#'
print.rapid_result <- function(x, type = "summary", ...) {
  if (type == "summary") {
    cat("RAPID Assessment\n")
    cat("================\n")
    cat("Method:", x$risk$method, "\n")
    cat("Risk level:", round(x$risk$confidence_rate * 100, 1), "%\n")
    cat("Records at risk:", x$risk$n_at_risk, "/",
        nrow(x$risk$rows_risk_df), "\n")

    # Threshold - check if categorical or numeric
    if (!is.null(x$risk$rows_risk_df$cat_tau)) {
      # Categorical
      cat("Threshold (tau):", x$risk$rows_risk_df$cat_tau[1], "\n")

    } else if (!is.null(x$risk$rows_risk_df$num_epsilon)) {
      # Numeric
      epsilon_val <- x$risk$rows_risk_df$num_epsilon[1]
      epsilon_scale <- x$risk$epsilon_scale

      if (!is.null(epsilon_scale) && !is.na(epsilon_scale)) {
        if (epsilon_scale == "percent") {
          cat("Threshold (epsilon):", epsilon_val * 100, "%\n")
        } else if (epsilon_scale == "proportion") {
          cat("Threshold (epsilon):", epsilon_val, "(",
              round(epsilon_val * 100, 1), "%)\n")
        }
      } else {
        # Absolute error or missing epsilon_scale
        cat("Threshold (epsilon):", epsilon_val, "\n")
      }
    }

  } else if (type == "high_risk") {
    at_risk <- x$risk$rows_risk_df[x$risk$rows_risk_df$at_risk == TRUE, , drop = FALSE]
    print(head(at_risk, 10))
    if (nrow(at_risk) > 10) {
      cat("... and", nrow(at_risk) - 10, "more\n")
    }
  }
  invisible(x)
}


#' ################################################################################
#' Plot method for RAPID results
#'
#' @param x A rapid_result object
#' @param tau_range Range of tau values for sensitivity curve (categorical)
#' @param epsilon_range Range of epsilon values for sensitivity curve (numeric)
#' @param ... Additional arguments (unused)
#' @export
plot.rapid_result <- function(x,
                              tau_range = seq(0, 1, by = 0.05),
                              epsilon_range = seq(0, 1, by = 0.05),
                              ...) {

  # Determine if categorical or numeric
  is_categorical <- x$risk$method %in% c("RCS_marginal", "RCS_conditional", "NCE")

  if (is_categorical) {
    # Categorical: RAPID vs tau
    if (x$risk$method == "RCS_marginal") {
      gains <- x$risk$normalized_gain
    } else if (x$risk$method == "RCS_conditional") {
      gains <- x$risk$relative_score
    } else if (x$risk$method == "NCE") {
      gains <- 1 - x$risk$normalized_ce  # Higher = more at risk
    }

    # Compare RAPID for different tau values
    rapid_values <- sapply(tau_range, function(tau) {
      sum(gains > tau) / length(gains)
    })

    # Plot
    plot(tau_range, rapid_values,
         type = "l", lwd = 2, col = "blue",
         xlab = expression(paste("Normalized Gain Threshold ", tau)),
         ylab = "RAPID (Proportion at Risk)",
         main = "Threshold Sensitivity Curve",
         ylim = c(0, 1))

    # Add current threshold
    current_tau <- x$risk$rows_risk_df$cat_tau[1]
    current_rapid <- x$risk$confidence_rate
    points(current_tau, current_rapid, pch = 19, col = "red", cex = 1.5)
    text(current_tau, current_rapid,
         paste0("tau = ", current_tau, "\nRAPID = ", round(current_rapid, 3)),
         pos = 4, col = "red")
    grid()

  } else {
    # Numeric: RAPID vs epsilon
    error_col <- grep("error", names(x$risk$rows_risk_df), value = TRUE)[1]
    errors <- x$risk$rows_risk_df[[error_col]]
    epsilon_scale <- x$risk$epsilon_scale

    # Compute RAPID for different epsilon values
    # (Re-evaluate existing errors with different thresholds to show sensitivity)
    if (!is.null(epsilon_scale) && !is.na(epsilon_scale) && epsilon_scale == "percent") {
      if (missing(epsilon_range)) {
        epsilon_range <- seq(0, 100, by = 5)
      }
      epsilon_range_prop <- epsilon_range / 100
      rapid_values <- sapply(epsilon_range_prop, function(eps) {
        sum(errors < eps) / length(errors)
      })

      # Plot with percent labels
      plot(epsilon_range, rapid_values,
           type = "l", lwd = 2, col = "blue",
           xlab = expression(paste("Error Threshold ", epsilon, " (%)")),
           ylab = "RAPID (Proportion at Risk)",
           main = "Threshold Sensitivity Curve",
           ylim = c(0, 1))

      # Current threshold in percent
      current_eps <- x$risk$rows_risk_df$num_epsilon[1] * 100

    } else {
      # Proportion or absolute scale
      rapid_values <- sapply(epsilon_range, function(eps) {
        sum(errors < eps) / length(errors)
      })

      # Determine label based on scale
      if (!is.null(epsilon_scale) && !is.na(epsilon_scale) && epsilon_scale == "proportion") {
        xlab_text <- expression(paste("Error Threshold ", epsilon, " (proportion)"))
      } else {
        # Absolute error - use actual unit
        xlab_text <- expression(paste("Error Threshold ", epsilon))
      }

      plot(epsilon_range, rapid_values,
           type = "l", lwd = 2, col = "blue",
           xlab = xlab_text,
           ylab = "RAPID (Proportion at Risk)",
           main = "Threshold Sensitivity Curve",
           ylim = c(0, 1))

      current_eps <- x$risk$rows_risk_df$num_epsilon[1]
    }

    current_rapid <- x$risk$confidence_rate
    points(current_eps, current_rapid, pch = 19, col = "red", cex = 1.5)
    text(current_eps, current_rapid,
         paste0("epsilon = ", round(current_eps, 3), "\nRAPID = ", round(current_rapid, 3)),
         pos = 4, col = "red")
    grid()
  }

  invisible(x)
}


#' ################################################################################
#' Summary method for RAPID results
#'
#' @param object A rapid_result object
#' @param ... Additional arguments (unused)
#' @export
summary.rapid_result <- function(object, ...) {
  cat("RAPID Risk Assessment Summary\n")
  cat("==============================\n\n")

  cat("Evaluation method:", object$risk$method, "\n")
  cat("Attacker model:", "Random Forest", "\n\n")  # Could be extracted from object

  cat("Risk Metrics:\n")
  cat("  Confidence rate:", round(object$risk$confidence_rate, 3), "\n")
  cat("  Records at risk:", object$risk$n_at_risk,
      "(", round(object$risk$percentage, 1), "%)\n\n")

  # Threshold info
  if (!is.null(object$risk$rows_risk_df$cat_tau)) {
    cat("  Threshold (tau):", object$risk$rows_risk_df$cat_tau[1], "\n\n")
  } else if (!is.null(object$risk$rows_risk_df$num_epsilon)) {
    cat("  Threshold (epsilon):", object$risk$rows_risk_df$num_epsilon[1], "\n\n")
  }

  cat("Model Performance:\n")
  if (!is.null(object$metrics$accuracy)) {
    cat("  Accuracy:", round(object$metrics$accuracy, 3), "\n")
  }
  if (!is.null(object$metrics$f1_score)) {
    cat("  F1-score:", round(object$metrics$f1_score, 3), "\n")
  }

  invisible(object)
}

#' #############################################################################
#' Print method for RAPID cross-validation results
#'
#' @param x A rapid_cv_result object
#' @param ... Additional arguments (unused)
print.rapid_cv_result <- function(x, ...) {
  cat("RAPID Cross-Validation Results\n")
  cat("===============================\n")
  cat("Evaluation method:", x$eval_method, "\n")
  cat("Attacker model:", x$model_type, "\n")
  cat("K-folds:", x$settings$k, "\n")

  # Always show both tau and epsilon
  if (x$settings$is_categorical) {
    cat("Threshold (tau):", ifelse(!is.na(x$threshold), x$threshold, "NA"), "\n")
    cat("Threshold (epsilon): NA\n")
  } else {
    cat("Threshold (tau): NA\n")
    cat("Threshold (epsilon):", ifelse(!is.na(x$threshold), x$threshold, "NA"), "\n")
  }

  cat("\nRisk Estimate:\n")
  cat("  Mean:", round(x$cv_summary$mean, 3), "\n")
  cat("  SD:", round(x$cv_summary$sd, 3), "\n")
  cat("  95% CI: [",
      round(x$cv_summary$ci_lower, 3), ", ",
      round(x$cv_summary$ci_upper, 3), "]\n", sep = "")

  invisible(x)
}
