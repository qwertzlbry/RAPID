#' @importFrom graphics plot
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
                              epsilon_range = seq(0, 20, by = 1),
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
         paste0("τ = ", current_tau, "\nRAPID = ", round(current_rapid, 3)),
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
         paste0("ε = ", round(current_eps, 3), "\nRAPID = ", round(current_rapid, 3)),
         pos = 4, col = "red")
    grid()
  }

  invisible(x)
}
