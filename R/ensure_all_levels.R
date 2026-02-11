#' @title Ensure Factor Levels Are Preserved for Model Fitting
#'
#' @description
#' This helper function ensures that all levels of the sensitive attribute in the original data—
#' including `"missing"` if it exists—are also present in the synthetic data. This prevents
#' model training functions (e.g., random forest or logistic regression) from silently dropping
#' levels that are absent in the synthetic dataset but exist in the original, particularly when
#' dealing with factor variables.
#'
#' If the `"missing"` level exists in \code{original_data[[sensitive_attribute]]} but is not present in
#' \code{synth[[sensitive_attribute]]}, the function adds a small number of dummy rows with that level
#' to the synthetic data to ensure consistency.
#'
#' @param synthetic_data A synthetic dataset (data frame).
#' @param original_data The original dataset (data frame).
#' @param sensitive_attribute The name (character string) of the sensitive variable (must be a factor).
#'
#' @return A modified version of \code{synthetic_data} with dummy rows added if necessary to ensure
#'         the `"missing"` level is preserved.
#'
#' @examples
#' set.seed(42)
#' original_data <- data.frame(sensitive = factor(c("A", "B", "missing")),
#'                     x = 1:3)
#' synthetic_data <- data.frame(sensitive = factor(c("A", "B")),
#'                     x = 10:11)
#' synthetic_data_new <- ensure_all_levels_used(synthetic_data, original_data, "sensitive")
#'
#' levels(synthetic_data_new$sensitive)  # includes "missing"
#' table(synthetic_data_new$sensitive)
#'
#' @export
ensure_all_levels_used <- function(synthetic_data, original_data, sensitive_attribute) {
  # Only proceed if sensitive_attribute is a factor and has a "missing" level in original_data
  if (is.factor(original_data[[sensitive_attribute]]) &&
      "missing" %in% levels(original_data[[sensitive_attribute]])) {

    # If "missing" level not present in synthetic_data, add dummy rows to inject it
    if (!"missing" %in% as.character(synthetic_data[[sensitive_attribute]])) {
      dummy_number <- max(1, ceiling(nrow(synthetic_data) * 0.001))  # 0.1% of rows, min 1
      dummy <- synthetic_data[rep(1, dummy_number), , drop = FALSE]
      dummy[[sensitive_attribute]] <- factor("missing", levels = levels(synthetic_data[[sensitive_attribute]]))

      synthetic_data <- rbind(synthetic_data, dummy)

      warning(sprintf(
        "Added %d dummy row(s) (%.2f%% of synthetic data) to ensure level 'missing' is present in '%s'.",
        dummy_number,
        100 * dummy_number / nrow(synthetic_data),
        sensitive_attribute
      ))

    } else {
      warning(sprintf(
        "No action taken: level 'missing' already present in synthetic variable '%s'.",
        sensitive_attribute
      ))
    }

  } else {
    warning(sprintf(
      "No action taken: variable '%s' is not a factor or does not include level 'missing' in original data.",
      sensitive_attribute
    ))
  }

  return(synthetic_data)
}
