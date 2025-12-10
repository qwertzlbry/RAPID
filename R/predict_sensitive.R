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

  predictions <- switch(model,

                        lm = {
                          stats::predict(fit, newdata = truth)
                        },

                        rf = {
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

                          predict(fit, data = truth, type = "response")$predictions
                        },

                        cart = {
                          if (!inherits(fit, "rpart")) {
                            stop("Model object is not of class 'rpart'")
                          }

                          is_classification <- is.factor(truth[[sensitive_var]])

                          if (is_classification) {
                            probs <- predict(fit, newdata = truth, type = "prob")
                            probs_df <- as.data.frame(probs)

                            # Get levels as CHARACTER
                            levels_response <- as.character(levels(truth[[sensitive_var]]))
                            colnames_probs <- colnames(probs_df)

                            # Check if all levels present
                            if (!all(levels_response %in% colnames_probs)) {
                              missing <- setdiff(levels_response, colnames_probs)
                              stop(sprintf("CART predictions missing levels: %s", paste(missing, collapse = ", ")))
                            }

                            # Return in correct order
                            probs_df[, levels_response, drop = FALSE]

                          } else {
                            predict(fit, newdata = truth)
                          }
                        },
                        gbm = {
                          if (!inherits(fit, "xgb.Booster")) {
                            stop("Model object is not of class 'xgb.Booster'")
                          }

                          # Reconstruct model matrix
                          X_test <- model.matrix(fit$.__x_formula__, data = truth)[, -1, drop = FALSE]
                          dtest  <- xgboost::xgb.DMatrix(X_test)

                          preds <- predict(fit, newdata = dtest)

                          y_true <- truth[[sensitive_var]]

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

                          probs <- predict(fit, newdata = truth, type = "response")

                          levels_response <- as.character(levels(truth[[sensitive_var]]))

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
