#' @title Evaluate Continuous Attribute Inference Risk (RAPID)
#'
#' @description
#' Computes RAPID disclosure risk for a continuous sensitive attribute by evaluating
#' how closely synthetic-trained models can predict true values. Risk is measured
#' by the proportion of records for which prediction errors fall below a chosen
#' absolute or relative threshold.
#'
#' @param A Vector of true sensitive values (numeric).
#' @param B Vector of predicted values (numeric), from a model trained on synthetic data.
#' @param original_data Data frame of original data (for reporting records with high-risk predictions).
#' @param num_error_metric Error metric: one of \code{"mae"}, \code{"rmse"}, \code{"rmae"}, or \code{"rrmse"}.
#' @param num_epsilon Numeric threshold (absolute error or percentage).
#' @param num_epsilon_type Threshold type: either \code{"Value"} for absolute error, or \code{"Percentage"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{rows_risk_n}{Number of records with small prediction error (at-risk cases).}
#'   \item{rows_risk_p}{Percentage of at-risk records relative to \code{original_data}.}
#'   \item{rows_risk_df}{Data frame of those records, including predicted and observed values and error metric.}
#' }
#'
#' @export
evaluate_numeric <- function(A, B, original_data,
                             num_error_metric,
                             num_epsilon,
                             num_epsilon_type,
                             return_all_records = FALSE) {

  num_epsilon_type <- match.arg(num_epsilon_type, choices = c("Value", "Percentage"))

  if (!is.numeric(num_epsilon)) {
    stop("`epsilon` must be numeric and represent either a raw error or percentage.")
  }

  # Compute absolute error-based metric
  measure <- switch(num_error_metric,
                    mae   = abs(A - B),
                    rmse  = sqrt((A - B)^2),
                    rmae  = abs(A - B) / mean(abs(A)),
                    rrmse = sqrt((A - B)^2) / stats::sd(A),
                    stop("Unsupported method. Use 'mae', 'rmse', 'rmae', or 'rrmse'.")
  )
  # Compute relative error for ALL records (for Percentage type)
  relative_error <- ifelse(
    A == 0, abs(B - A) / (abs(A) + 1e-6) * 100,
    abs(B - A) / abs(A) * 100
  )

  # Evaluate threshold match
  switch(num_epsilon_type,

         Value = {
           at_risk <- measure < num_epsilon
           metric_col <- measure
           metric_name <- num_error_metric
         },

         Percentage = {
           # Relative error: avoid division by zero

           at_risk <- relative_error < num_epsilon
           metric_col <- relative_error
           metric_name <- "relative_error_percent"
         }
  )

  # Output reporting table
  result <- data.frame(
    original_data,
    Original = A,
    Predicted = B,
    error_metric = metric_col,
    num_epsilon = num_epsilon,
    at_risk = at_risk,
    stringsAsFactors = FALSE
  )
  colnames(result)[colnames(result) == "error_metric"] <- metric_name

  # Filter to at-risk only if needed
  if (return_all_records) {
    result_df <- result
  } else {
    result_df <- result[result$at_risk == TRUE, ]
  }

  return(list(
    confidence_rate = sum(at_risk)/ nrow(original_data),
    n_at_risk = sum(at_risk),
    percentage = 100 * sum(at_risk) / nrow(original_data),
    rows_risk_df = result_df
  ))
}
