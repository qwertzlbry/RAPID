#' @title Fit Predictive Model on Synthetic Data
#'
#' @description
#' Trains a regression or classification model on synthetic data as part of the RAPID disclosure risk workflow.
#' Supports linear models (\code{"lm"}), random forests via \pkg{ranger}, or CART via \pkg{rpart}.
#' Automatically enables probability output for classification when the sensitive attribute is categorical.
#'
#' @param model A string specifying the model type: \code{"lm"}, \code{"rf"}, or \code{"cart"}.
#' @param formula A model formula of the form \code{y ~ x1 + x2 + ...}.
#' @param synthetic_data A \code{data.frame} containing the synthetic data for training.
#' @param original_data A \code{data.frame} used to check whether the sensitive attribute is categorical (for probability output).
#' @param sensitive_attribute A string; name of the sensitive attribute (outcome) to be predicted.
#' @param lm.control Optional list of arguments passed to \code{lm()}.
#' @param ranger.control Optional list of arguments passed to \code{ranger()}.
#' @param rpart.control Optional list of arguments passed to \code{rpart()}.
#'
#' @return A fitted model object from \code{lm()}, \code{ranger()}, or \code{rpart()}.
#' @export
fit_model <- function(model, formula, synthetic_data, original_data, sensitive_attribute,
                      lm.control = list(),
                      ranger.control = list(),
                      rpart.control = list(),
                      xgb.control = list()) {

  # Auto-convert character to factor for classification
  if (is.character(synthetic_data[[sensitive_attribute]])) {
    synthetic_data[[sensitive_attribute]] <- factor(synthetic_data[[sensitive_attribute]])
  }
  if (is.character(original_data[[sensitive_attribute]])) {
    original_data[[sensitive_attribute]] <- factor(original_data[[sensitive_attribute]])
  }


  fit <- switch(model,

                lm = {
                  args <- c(list(formula = formula, data = synthetic_data), lm.control)
                  do.call(stats::lm, args)
                },

                rf = {
                  use_prob <- is.factor(original_data[[sensitive_attribute]])
                  args <- c(
                    list(
                      formula = formula,
                      data = synthetic_data,
                      probability = use_prob
                    ),
                    ranger.control
                  )
                  do.call(ranger::ranger, args)
                },

                cart = {
                  # Determine method based on response type
                  method_type <- if (is.factor(synthetic_data[[sensitive_attribute]])) {
                    "class"
                  } else {
                    "anova"
                  }

                  args <- c(
                    list(
                      formula = formula,
                      data = synthetic_data,
                      method = method_type
                    ),
                    rpart.control
                  )
                  do.call(rpart::rpart, args)
                },

                gbm = {

                  y_raw <- synthetic_data[[sensitive_attribute]]

                  # DESIGN MATRIX

                  X <- model.matrix(formula, data = synthetic_data)[, -1, drop = FALSE]

                  # REGRESSION
                  if (!is.factor(y_raw)) {

                    dtrain <- xgboost::xgb.DMatrix(
                      data  = X,
                      label = y_raw
                    )

                    params <- c(
                      list(
                        objective = "reg:squarederror",
                        eval_metric = "rmse",
                        verbosity = 0
                      ),
                      xgb.control
                    )

                    model <- xgboost::xgb.train(
                      params  = params,
                      data    = dtrain,
                      nrounds = 100,
                      verbose = 0
                    )

                    model$.__x_formula__ <- formula
                    return(model)
                  }

                  # CLASSIFICATION
                  y_factor  <- factor(y_raw)
                  y_encoded <- as.integer(y_factor) - 1
                  K <- length(levels(y_factor))

                  dtrain <- xgboost::xgb.DMatrix(
                    data  = X,
                    label = y_encoded
                  )

                  # ----- BINARY -----
                  if (K == 2) {

                    params <- c(
                      list(
                        objective = "binary:logistic",
                        eval_metric = "logloss",
                        verbosity = 0
                      ),
                      xgb.control
                    )

                    model <- xgboost::xgb.train(
                      params  = params,
                      data    = dtrain,
                      nrounds = 100,
                      verbose = 0
                    )

                    # ----- MULTICLASS -----
                  } else {

                    params <- c(
                      list(
                        objective = "multi:softprob",
                        num_class = K,
                        eval_metric = "mlogloss",
                        verbosity = 0
                      ),
                      xgb.control
                    )

                    model <- xgboost::xgb.train(
                      params  = params,
                      data    = dtrain,
                      nrounds = 100,
                      verbose = 0
                    )
                  }

                  model$.__x_formula__ <- formula
                  return(model)
                },
                logit = {
                  if (!is.factor(synthetic_data[[sensitive_attribute]])) {
                    stop("Logit requires a categorical sensitive attribute.")
                  }

                  glm(formula = formula, data = synthetic_data, family = binomial())
                },
                stop("Unsupported model type. Use 'lm', 'rf', or 'cart'.")
  )

  return(fit)
}
