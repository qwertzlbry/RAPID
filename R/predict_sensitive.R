#' Predict Sensitive Attribute from Fitted Model
#'
#' @param model_type A string: either "lm" or "rf".
#' @param fit The fitted model object.
#' @param original_data The original data (data frame) on which the predictions are made.
#' @param sensitive_attribute The name of the sensitive attribute (string).
#'
#' @return `B` — predicted values (vector for continuous, matrix or factor for categorical).
#' @export
predict_sensitive <- function(model, fit, truth, sensitive_var) {
  if (model == "lm") {
    return(stats::predict(fit, newdata = truth))
  }
  
  if (model == "rf") {
    if (!inherits(fit, "ranger")) {
      stop("Model object is not of class 'ranger'")
    }
    # Ensure categorical target remains a factor
    if (is.factor(fit$call$data[[as.character(fit$call$formula[[2]])]])) {
      # Classification: force truth var to be factor (if needed)
      if (!is.factor(truth[[sensitive_var]])) {
        truth[[sensitive_var]] <- factor(truth[[sensitive_var]],
                                         levels = levels(fit$forest$independent.variable.names))
      }
    }
    
    return(predict(fit, data = truth, type = "response")$predictions)
  }
  
  stop("Unsupported model type.")
}