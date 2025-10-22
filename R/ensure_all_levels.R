#' @title Ensure Factor Levels Are Preserved for Model Fitting
#'
#' @description
#' This helper function ensures that all levels of the sensitive attribute in the original data—
#' including `"missing"` if it exists—are also present in the synthetic data. This prevents
#' model training functions (e.g., random forest or logistic regression) from silently dropping
#' levels that are absent in the synthetic dataset but exist in the original, particularly when
#' dealing with factor variables.
#'
#' If the `"missing"` level exists in \code{truth[[sensitive_var]]} but is not present in
#' \code{synth[[sensitive_var]]}, the function adds a small number of dummy rows with that level
#' to the synthetic data to ensure consistency.
#'
#' @param synth A synthetic dataset (data frame).
#' @param truth The original dataset (data frame).
#' @param sensitive_var The name (character string) of the sensitive variable (must be a factor).
#'
#' @return A modified version of \code{synth} with dummy rows added if necessary to ensure
#'         the `"missing"` level is preserved.
#'
#' @examples
#' set.seed(42)
#' truth <- data.frame(sensitive = factor(c("A", "B", "missing")),
#'                     x = 1:3)
#' synth <- data.frame(sensitive = factor(c("A", "B")),
#'                     x = 10:11)
#' synth_new <- ensure_all_levels_used(synth, truth, "sensitive")
#'
#' levels(synth_new$sensitive)  # includes "missing"
#' table(synth_new$sensitive)
#'
#' @export
ensure_all_levels_used <- function(synth, truth, sensitive_var) {
  # Only proceed if sensitive_var is a factor and has a "missing" level in truth
  if (is.factor(truth[[sensitive_var]]) &&
      "missing" %in% levels(truth[[sensitive_var]])) {
    
    # If "missing" level not present in synth, add dummy rows to inject it
    if (!"missing" %in% as.character(synth[[sensitive_var]])) {
      dummy_number <- max(1, ceiling(nrow(synth) * 0.001))  # 0.1% of rows, min 1
      dummy <- synth[rep(1, dummy_number), , drop = FALSE]
      dummy[[sensitive_var]] <- factor("missing", levels = levels(synth[[sensitive_var]]))
      
      synth <- rbind(synth, dummy)
      
      warning(sprintf(
        "Added %d dummy row(s) (%.2f%% of synthetic data) to ensure level 'missing' is present in '%s'.",
        dummy_number,
        100 * dummy_number / nrow(synth),
        sensitive_var
      ))
      
    } else {
      warning(sprintf(
        "No action taken: level 'missing' already present in synthetic variable '%s'.",
        sensitive_var
      ))
    }
    
  } else {
    warning(sprintf(
      "No action taken: variable '%s' is not a factor or does not include level 'missing' in original data.",
      sensitive_var
    ))
  }
  
  return(synth)
}