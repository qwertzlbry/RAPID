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
#' @param truth Data frame of original data (for reporting records with high-risk predictions).
#' @param method Error metric: one of \code{"mae"}, \code{"rmse"}, \code{"rmae"}, or \code{"rrmse"}.
#' @param epsilon Numeric threshold (absolute error or percentage).
#' @param epsilon_type Threshold type: either \code{"Value"} for absolute error, or \code{"Percentage"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{rows_risk_n}{Number of records with small prediction error (at-risk cases).}
#'   \item{rows_risk_p}{Percentage of at-risk records relative to \code{truth}.}
#'   \item{rows_risk_df}{Data frame of those records, including predicted and observed values and error metric.}
#' }
#'
#' @export
evaluate_numeric <- function(A, B, truth,
                             method,
                             epsilon,
                             epsilon_type) {
  
  epsilon_type <- match.arg(epsilon_type, choices = c("Value", "Percentage"))
  
  if (!is.numeric(epsilon)) {
    stop("`epsilon` must be numeric and represent either a raw error or percentage.")
  }
  
  # Compute absolute error-based metric
  measure <- switch(method,
                    mae   = abs(A - B),
                    rmse  = sqrt((A - B)^2),
                    rmae  = abs(A - B) / mean(abs(A)),
                    rrmse = sqrt((A - B)^2) / stats::sd(A),
                    stop("Unsupported method. Use 'mae', 'rmse', 'rmae', or 'rrmse'.")
  )
  
  # Evaluate threshold match
  switch(epsilon_type,
         
         Value = {
           idx <- which(measure < epsilon)
           label <- method
           extra_col <- measure[idx]
         },
         
         Percentage = {
           # Relative error: avoid division by zero
           relative_error <- ifelse(
             A == 0, abs(B - A) / (abs(A) + 1e-6) * 100,
             abs(B - A) / abs(A) * 100
           )
           idx <- which(relative_error < epsilon)
           label <- "relative_error_percent"
           extra_col <- relative_error[idx]
         }
  )
  
  # Output reporting table
  result <- data.frame(truth[idx, , drop = FALSE],
                       Original = A[idx],
                       Predicted = B[idx],
                       Dist.Metric = extra_col)
  colnames(result)[ncol(result)] <- label
  
  return(list(
    rows_risk_n = length(idx),
    rows_risk_p = 100 * length(idx) / nrow(truth),
    rows_risk_df = result
  ))
}