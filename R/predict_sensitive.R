#' Predict Sensitive Attribute from Fitted Model
#'
#' @param model_type A string: either "lm" or "rf".
#' @param fit The fitted model_type object.
#' @param original_data The original data (data frame) on which the predictions are made.
#' @param sensitive_attribute The name of the sensitive attribute (string).
#'
#' @return `B` — predicted values (vector for continuous, matrix or factor for categorical).
#' @export
predict_sensitive <- function(model_type, fit, original_data, sensitive_attribute) {

  predictions <- switch(model_type,

                        lm = {
                          stats::predict(fit, newdata = original_data)
                        },

                        rf = {
                          if (!inherits(fit, "ranger")) {
                            stop("Model object is not of class 'ranger'")
                          }

                          # Ensure categorical target remains a factor
                          if (is.factor(fit$call$data[[as.character(fit$call$formula[[2]])]])) {
                            # Classification: force original_data var to be factor (if needed)
                            if (!is.factor(original_data[[sensitive_attribute]])) {
                              original_data[[sensitive_attribute]] <- factor(original_data[[sensitive_attribute]],
                                                               levels = levels(fit$forest$independent.variable.names))
                            }
                          }

                          predict(fit, data = original_data, type = "response")$predictions
                        },

                        cart = {
                          if (!inherits(fit, "rpart")) {
                            stop("Model object is not of class 'rpart'")
                          }

                          is_classification <- is.factor(original_data[[sensitive_attribute]])

                          if (is_classification) {
                            probs <- predict(fit, newdata = original_data, type = "prob")
                            probs_df <- as.data.frame(probs)

                            # Get levels as CHARACTER
                            levels_response <- as.character(levels(original_data[[sensitive_attribute]]))
                            colnames_probs <- colnames(probs_df)

                            # Check if all levels present
                            if (!all(levels_response %in% colnames_probs)) {
                              missing <- setdiff(levels_response, colnames_probs)
                              stop(sprintf("CART predictions missing levels: %s", paste(missing, collapse = ", ")))
                            }

                            # Return in correct order
                            probs_df[, levels_response, drop = FALSE]

                          } else {
                            predict(fit, newdata = original_data)
                          }
                        },
                        gbm = {
                          if (!inherits(fit, "xgb.Booster")) {
                            stop("Model object is not of class 'xgb.Booster'")
                          }

                          # Reconstruct model matrix
                          X_test <- model.matrix(fit$.__x_formula__, data = original_data)[, -1, drop = FALSE]
                          dtest  <- xgboost::xgb.DMatrix(X_test)

                          preds <- predict(fit, newdata = dtest)

                          y_true <- original_data[[sensitive_attribute]]

                          ## REGRESSION
                          if (!is.factor(y_true)) {
                            return(preds)
                          }

                          ## CLASSIFICATION

                          levels_response <- as.character(levels(y_true))
                          K <- length(levels_response)

                          ## ----- BINARY -----
                          if (K == 2) {
                            probs_df <- data.frame(
                              class1 = 1 - preds,
                              class2 = preds
                            )
                            colnames(probs_df) <- levels_response
                            return(probs_df)
                          }

                          ## ----- MULTICLASS -----
                          probs_mat <- matrix(preds, ncol = K, byrow = TRUE)
                          probs_df  <- as.data.frame(probs_mat)
                          colnames(probs_df) <- levels_response
                          return(probs_df)
                        },

                        logit = {
                          if (!inherits(fit, "glm")) {
                            stop("Model object is not of class 'glm'")
                          }

                          probs <- predict(fit, newdata = original_data, type = "response")

                          levels_response <- as.character(levels(original_data[[sensitive_attribute]]))

                          if (length(levels_response) != 2) {
                            stop("Logit only supports binary classification.")
                          }

                          probs_df <- data.frame(
                            class1 = 1 - probs,
                            class2 = probs
                          )
                          colnames(probs_df) <- levels_response
                          probs_df
                        },


                        stop("Unsupported model type. Use 'lm', 'rf', or 'cart'.")
  )

  return(predictions)
}
