#' Handle Missing Values in the Sensitive Attribute
#'
#' @description
#' Deals with missing values in the sensitive attribute for both original and synthetic data.
#' - For factors: replaces NAs with a new level "missing".
#' - For numerics: allows "drop", "constant", or "median" strategies.
#'
#' @param original_data The original dataset (data frame).
#' @param synthetic_data The synthetic dataset (data frame).
#' @param sensitive_attribute Name of the sensitive attribute (string).
#' @param num_na_strategy Strategy for handling missing numeric values: one of "drop", "constant", or "median".
#' @param num_constant_value Value used for imputation if \code{num_na_strategy = "constant"}.
#'
#' @return A list with cleaned \code{original_data} and \code{synthetic_data}.
#' @export
handle_sensitive_na <- function(original_data,
                                synthetic_data,
                                sensitive_attribute,
                                num_na_strategy = c("constant", "drop", "median"),
                                num_constant_value = 0) {

  num_na_strategy <- match.arg(num_na_strategy)

  # Handle factor (categorical) sensitive attribute
  if (is.factor(original_data[[sensitive_attribute]]) &&
      (anyNA(original_data[[sensitive_attribute]]) || anyNA(synthetic_data[[sensitive_attribute]]))) {

    original_data[[sensitive_attribute]] <- addNA(original_data[[sensitive_attribute]])
    synthetic_data[[sensitive_attribute]] <- addNA(synthetic_data[[sensitive_attribute]])

    levels(original_data[[sensitive_attribute]])[is.na(levels(original_data[[sensitive_attribute]]))] <- "missing"
    levels(synthetic_data[[sensitive_attribute]])[is.na(levels(synthetic_data[[sensitive_attribute]]))] <- "missing"

    warning(sprintf("Missing values in categorical sensitive attribute '%s': replaced with 'missing' level.", sensitive_attribute))

    # Handle numeric sensitive attribute
  } else if (is.numeric(original_data[[sensitive_attribute]]) &&
             (anyNA(original_data[[sensitive_attribute]]) || anyNA(synthetic_data[[sensitive_attribute]]))) {

    if (num_na_strategy == "drop") {
      warning(sprintf("Missing values in numeric sensitive attribute '%s': rows with NA dropped.", sensitive_attribute))
      keep <- complete.cases(original_data[[sensitive_attribute]], synthetic_data[[sensitive_attribute]])
      original_data <- original_data[keep, , drop = FALSE]
      synthetic_data <- synthetic_data[keep, , drop = FALSE]

    } else if (num_na_strategy == "constant") {
      warning(sprintf("Missing values in numeric sensitive attribute '%s': replaced with constant value = %s.",
                      sensitive_attribute, as.character(num_constant_value)))
      original_data[[sensitive_attribute]][is.na(original_data[[sensitive_attribute]])] <- num_constant_value
      synthetic_data[[sensitive_attribute]][is.na(synthetic_data[[sensitive_attribute]])] <- num_constant_value

    } else if (num_na_strategy == "median") {
      warning(sprintf("Missing values in numeric sensitive attribute '%s': replaced with median from each dataset.", sensitive_attribute))
      orig_med <- median(original_data[[sensitive_attribute]], na.rm = TRUE)
      synth_med <- median(synthetic_data[[sensitive_attribute]], na.rm = TRUE)
      original_data[[sensitive_attribute]][is.na(original_data[[sensitive_attribute]])] <- orig_med
      synthetic_data[[sensitive_attribute]][is.na(synthetic_data[[sensitive_attribute]])] <- synth_med
    }
  }

  return(list(
    original_data = original_data,
    synthetic_data = synthetic_data
  ))
}
