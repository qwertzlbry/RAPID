#' @title Evaluate Continuous Attribute Inference Risk (RAPID)
#'
#' @description
#' Computes RAPID disclosure risk for a continuous sensitive attribute by evaluating
#' how closely synthetic-trained models can predict true values. Risk is measured
#' by the proportion of records for which prediction errors fall below a chosen
#' threshold.
#'
#' @param A Vector of true sensitive values (numeric).
#' @param B Vector of predicted values (numeric), from a model trained on synthetic data.
#' @param original_data Data frame of original data (for reporting records with high-risk predictions).
#' @param num_error_metric Error metric: one of \code{"symmetric"} (default, recommended),
#'   \code{"stabilised_relative"}, or \code{"absolute"}. See Details.
#' @param num_epsilon Numeric threshold. Interpretation depends on \code{num_error_metric}.
#'   For relative metrics, this is a relative error threshold; for \code{"absolute"}, this is
#'   an absolute error threshold on the original scale.
#' @param num_epsilon_scale Scale of \code{num_epsilon} for relative metrics:
#'   \code{"proportion"} means \code{0.05} corresponds to 5\%, and \code{"percent"} means \code{5}
#'   corresponds to 5\%. Only used when \code{num_error_metric} is \code{"symmetric"} or
#'   \code{"stabilised_relative"}.
#' @param num_delta Smoothing constant to prevent division by zero for relative metrics
#'   (default 0.01). Only used for \code{"symmetric"} and \code{"stabilised_relative"}.
#' @param return_all_records Logical; if \code{TRUE}, return all records with their risk
#'   status. If \code{FALSE} (default), return only at-risk records.
#'
#' @details
#' Three error metrics are available:
#' \itemize{
#'   \item \code{symmetric}: Symmetric percentage error (recommended). Treats predicted and
#'     true values symmetrically:
#'     \eqn{2|y - \hat{y}| / (|y| + |\hat{y}| + 2\delta)}.
#'   \item \code{stabilised_relative}: Stabilised relative error:
#'     \eqn{|y - \hat{y}| / (|y| + \delta)}.
#'   \item \code{absolute}: Absolute error:
#'     \eqn{|y - \hat{y}|}.
#' }
#'
#' For relative metrics (\code{"symmetric"} and \code{"stabilised_relative"}), errors are computed
#' on a 0--1 scale. The threshold \code{num_epsilon} is interpreted according to
#' \code{num_epsilon_scale}:
#' \itemize{
#'   \item \code{num_epsilon_scale = "proportion"}: use \code{0.05} for 5\%
#'   \item \code{num_epsilon_scale = "percent"}: use \code{5} for 5\% (internally converted to 0.05)
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{confidence_rate}{Proportion of records at risk (between 0 and 1).}
#'   \item{n_at_risk}{Number of records with prediction errors below the threshold.}
#'   \item{percentage}{Percentage of at-risk records (0--100 scale).}
#'   \item{rows_risk_df}{Data frame of records (all or at-risk only, depending on
#'     \code{return_all_records}) with original values, predictions, error metrics,
#'     and risk status. For relative metrics the error column is
#'     \code{symmetric_error} or \code{stabilised_relative_error}; for absolute it is
#'     \code{absolute_error}. The table also includes \code{num_epsilon_input} and the
#'     converted \code{num_epsilon} actually used for the comparison.}
#' }
#'
#' @examples
#' \dontrun{
#' # Relative metric (symmetric), epsilon given in percent (5 means 5%)
#' result <- evaluate_numeric(
#'   A = original_data$income,
#'   B = predictions,
#'   original_data = original_data,
#'   num_error_metric = "symmetric",
#'   num_epsilon = 5,
#'   num_epsilon_scale = "percent"
#' )
#'
#' # Relative metric, epsilon given as proportion (0.05 means 5%)
#' result_prop <- evaluate_numeric(
#'   A = original_data$income,
#'   B = predictions,
#'   original_data = original_data,
#'   num_error_metric = "symmetric",
#'   num_epsilon = 0.05,
#'   num_epsilon_scale = "proportion"
#' )
#'
#' # Absolute threshold (e.g., income in dollars)
#' result_abs <- evaluate_numeric(
#'   A = original_data$income,
#'   B = predictions,
#'   original_data = original_data,
#'   num_error_metric = "absolute",
#'   num_epsilon = 1000,
#'   return_all_records = TRUE
#' )
#' }
#'
#' @export

evaluate_numeric <- function(A, B, original_data,
                             num_error_metric = c("symmetric", "stabilised_relative", "absolute"),
                             num_epsilon,
                             num_epsilon_scale = c("proportion", "percent"), # 0.05 vs 5
                             num_delta = 0.01,
                             return_all_records = FALSE) {

  num_error_metric <- match.arg(num_error_metric)
  num_epsilon_scale <- match.arg(num_epsilon_scale)


  # Validation
  if (!is.logical(return_all_records) || length(return_all_records) != 1) {
    stop("`return_all_records` must be a single TRUE/FALSE.")
  }

  if (!is.numeric(num_epsilon)) {
    stop("`num_epsilon` must be numeric and represent either a raw error or percentage.")
  }
  if (!is.numeric(num_delta) || num_delta <= 0) {
    stop("num_delta must be a positive numeric value")
  }
  if (!is.numeric(A) || !is.numeric(B)) stop("`A` and `B` must be numeric.")
  if (length(A) != length(B)) stop("`A` and `B` must have the same length.")
  if (!is.data.frame(original_data)) stop("`original_data` must be a data.frame.")
  if (nrow(original_data) != length(A)) stop("nrow(original_data) must equal length(A).")


   is_relative <- num_error_metric %in% c("symmetric","stabilised_relative")
   if (is_relative && num_epsilon_scale == "proportion" && num_epsilon > 1) {
     warning("`num_epsilon` is > 1 but `num_epsilon_scale='proportion'`. Did you mean `num_epsilon_scale='percent'` (e.g., 5 for 5%)?")
   }

   if (is_relative) {

  # Compute error based on chosen metric
  error_values <- switch(num_error_metric,
                         symmetric = {
                          2 * abs(A - B) / (abs(A) + abs(B) + 2*num_delta)
                         },
                         stabilised_relative = {
                         abs(A - B) / (abs(A) + num_delta)
                         }
  )
  metric_name <- paste0(num_error_metric, "_error")
  epsilon_used <- if (num_epsilon_scale == "percent") num_epsilon / 100 else num_epsilon
   } else {
    # num_epsilon_type == "absolute"
    # Absolute error on original scale
    error_values <- abs(A - B)
    epsilon_used <- num_epsilon
    metric_name <- "absolute_error"

  }

  at_risk <- error_values < epsilon_used

  # Output reporting table
  result <- data.frame(
    original_data,
    Original = A,
    Predicted = B,
    error_metric = error_values,
    num_epsilon_input = num_epsilon,
    num_epsilon = epsilon_used,
    num_epsilon_scale = if (is_relative) num_epsilon_scale else NA_character_,
    at_risk = at_risk,
    stringsAsFactors = FALSE
  )
  colnames(result)[colnames(result) == "error_metric"] <- metric_name

  # Filter to at-risk only if return_all_records is FALSE
  if (!return_all_records) {
    result <- result[result$at_risk == TRUE, , drop = FALSE]
  }

  return(list(
    method = num_error_metric,
    confidence_rate = sum(at_risk)/ nrow(original_data),
    n_at_risk = sum(at_risk),
    percentage = 100 * sum(at_risk) / nrow(original_data),
    rows_risk_df = result
  ))
}
