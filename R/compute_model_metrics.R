#' @title Compute Basic Model Evaluation Metrics
#'
#' @description
#' Computes evaluation metrics for classification (accuracy) or regression tasks 
#' (MAE, RMSE, relative MAE, relative RMSE), used for assessing attacker model performance 
#' in the RAPID disclosure risk framework.
#'
#' @param A True values from the original dataset (numeric for continuous, factor for categorical).
#' @param B Predicted values from the attacker model. Should be:
#'   - A numeric vector for continuous attributes.
#'   - A factor vector or probability matrix (with class labels as column names) for categorical attributes.
#'
#' @return A named list with one or more of the following elements:
#'   \item{accuracy}{Classification accuracy (for categorical variables).}
#'   \item{mae}{Mean Absolute Error (for continuous variables).}
#'   \item{rmse}{Root Mean Squared Error (for continuous variables).}
#'   \item{rmae}{Relative MAE: \code{mae / mean(abs(A))}.}
#'   \item{rrmse}{Relative RMSE: \code{rmse / sd(A)}.}
#'
#' @examples
#' A_num <- c(100, 200, 300)
#' B_num <- c(110, 190, 310)
#' compute_model_metrics(A_num, B_num)
#'
#' A_cat <- factor(c("yes", "no", "yes"))
#' B_mat <- matrix(c(0.7, 0.3,
#'                   0.4, 0.6,
#'                   0.8, 0.2), nrow = 3, byrow = TRUE)
#' colnames(B_mat) <- c("yes", "no")
#' compute_model_metrics(A_cat, B_mat)
#'
#' @export
compute_model_metrics <- function(A, B) {
  if (is.factor(A)) {
    # Classification case
    if (is.matrix(B)) {
      pred_class <- colnames(B)[apply(B, 1, which.max)]
    } else {
      pred_class <- B
    }
    accuracy <- mean(pred_class == A)
    return(list(accuracy = accuracy))
  } else {
    # Regression case
    mae <- mean(abs(A - B))
    rmse <- sqrt(mean((A - B)^2))
    rmae <- mae / mean(abs(A))
    rrmse <- rmse / sd(A)
    return(list(mae = mae, rmse = rmse, rmae = rmae, rrmse = rrmse))
  }
}