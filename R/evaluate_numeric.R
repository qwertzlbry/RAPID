#' @title Evaluate Continuous Attribute Inference Risk (RAPID)
#'
#' @description
#' Computes RAPID disclosure risk for a continuous sensitive attribute by evaluating
#' how closely synthetic-trained models can predict true values. Risk is measured
#' by the proportion of records for which prediction errors fall below a chosen
#' threshold. Two threshold types are supported: percentage-based (relative error)
#' and absolute (raw error on the original scale).
#'
#' @param A Vector of true sensitive values (numeric).
#' @param B Vector of predicted values (numeric), from a model trained on synthetic data.
#' @param original_data Data frame of original data (for reporting records with high-risk predictions).
#' @param num_error_metric Error metric: one of \code{"symmetric"} (default, recommended),
#'   \code{"stabilised_relative"}, or \code{"absolute"}. See Details.
#' @param num_epsilon Numeric threshold. For percentage-based metrics, specify as a
#'   percentage (e.g., 5 for 5\%). For absolute error, specify in the units of the
#'   sensitive attribute (e.g., 1000 for $1000).
#' @param num_epsilon_type Threshold type: either \code{"percentage"} for relative error
#'   metrics, or \code{"absolute"} for raw error on the original scale. Default is \code{"percentage"}.
#' @param delta Smoothing constant to prevent division by zero for percentage-based metrics
#'   (default 0.01). Only used when \code{num_epsilon_type = "percentage"}.
#' @param return_all_records Logical; if \code{TRUE}, return all records with their risk
#'   status. If \code{FALSE} (default), return only at-risk records.
#'
#' @details
#' Three error metrics are available:
#' \itemize{
#'   \item \code{symmetric}: Symmetric percentage error (recommended). Treats predicted and
#'     true values symmetrically: \eqn{2|y - \hat{y}| / (|y| + |\hat{y}| + 2\delta)}.
#'   \item \code{stabilised_relative}: Stabilised relative error: \eqn{|y - \hat{y}| / (|y| + \delta)}.
#'   \item \code{absolute}: Absolute error: \eqn{|y - \hat{y}|}. Only meaningful with
#'     \code{num_epsilon_type = "absolute"}.
#' }
#'
#' When \code{num_epsilon_type = "percentage"}, errors are expressed as percentages (0-100 scale).
#' When \code{num_epsilon_type = "absolute"}, raw absolute errors are used.
#'
#' @return A list containing:
#' \describe{
#'   \item{confidence_rate}{Proportion of records at risk (between 0 and 1).}
#'   \item{n_at_risk}{Number of records with prediction errors below the threshold.}
#'   \item{percentage}{Percentage of at-risk records.}
#'   \item{rows_risk_df}{Data frame of records (all or at-risk only, depending on
#'     \code{return_all_records}) with original values, predictions, error metrics,
#'     and risk status.}
#' }
#'
#' @examples
#' \dontrun{
#' # Percentage-based with symmetric error (recommended)
#' result <- evaluate_numeric(
#'   A = original_data$income,
#'   B = predictions,
#'   original_data = original_data,
#'   num_error_metric = "symmetric",
#'   num_epsilon = 5,  # 5% error threshold
#'   num_epsilon_type = "percentage"
#' )
#'
#' # Absolute threshold
#' result <- evaluate_numeric(
#'   A = original_data$income,
#'   B = predictions,
#'   original_data = original_data,
#'   num_epsilon = 1000,  # $1000 error threshold
#'   num_epsilon_type = "absolute"
#' )
#' }
#'
#' @export
evaluate_numeric <- function(A, B, original_data,
                             num_error_metric = c("symmetric", "stabilised_relative", "absolute"),
                             num_epsilon,
                             num_epsilon_type = c("percentage", "absolute"),
                             num_delta = 0.01,
                             error_as_percentage = FALSE,
                             return_all_records = FALSE) {

  num_error_metric <- match.arg(num_error_metric)

  if (!is.numeric(num_epsilon)) {
    stop("`epsilon` must be numeric and represent either a raw error or percentage.")
  }
  if (!is.numeric(num_delta) || num_delta <= 0) {
    stop("delta must be a positive numeric value")
  }


  if (num_epsilon_type == "percentage") {

  # Compute error based on chosen metric
  error_values <- switch(num_error_metric,
                         symmetric = {
                          err <-  2 * abs(A - B) / (abs(A) + abs(B) + 2*num_delta)
                          if (error_as_percentage) err * 100 else err
                         },
                         stabilised_relative = {
                          err <- abs(A - B) / (abs(A) + num_delta)
                           if (error_as_percentage) err * 100 else err
                         },
                         absolute = {
                           stop("Percentage-based threshold is not meaningful for absolute error metric. Use num_epsilon_type = 'absolute'.")
                         }

  )
  metric_name <- paste0(num_error_metric, "_error_pct")

  } else {
    # num_epsilon_type == "absolute"
    # Absolute error on original scale
    error_values <- abs(A - B)
    metric_name <- "absolute_error"

  }

  at_risk <- error_values < num_epsilon

  # Output reporting table
  result <- data.frame(
    original_data,
    Original = A,
    Predicted = B,
    error_metric = error_values,
    num_epsilon = num_epsilon,
    num_epsilon_type = num_epsilon_type,
    at_risk = at_risk,
    stringsAsFactors = FALSE
  )
  colnames(result)[colnames(result) == "error_metric"] <- metric_name

  # Filter to at-risk only if return_all_records is FALSE
  if (!return_all_records) {
    result <- result[result$at_risk == TRUE, ]
  }

  return(list(
    confidence_rate = sum(at_risk)/ nrow(original_data),
    n_at_risk = sum(at_risk),
    percentage = 100 * sum(at_risk) / nrow(original_data),
    rows_risk_df = result
  ))
}
