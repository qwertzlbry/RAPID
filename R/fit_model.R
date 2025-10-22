#' @title Fit Predictive Model on Synthetic Data
#'
#' @description
#' Trains a regression or classification model on synthetic data as part of the RAPID disclosure risk workflow.
#' Supports linear models (\code{"lm"}) or random forests via \pkg{ranger}. 
#' Automatically enables probability output for classification when the sensitive attribute is categorical.
#'
#' @param model A string specifying the model type: either \code{"lm"} or \code{"rf"}.
#' @param formula A model formula of the form \code{y ~ x1 + x2 + ...}.
#' @param synthetic_data A \code{data.frame} containing the synthetic data for training.
#' @param original_data A \code{data.frame} used to check whether the sensitive attribute is categorical (for probability output).
#' @param sensitive_attribute A string; name of the sensitive attribute (outcome) to be predicted.
#' @param lm.control Optional list of arguments passed to \code{lm()}.
#' @param ranger.control Optional list of arguments passed to \code{ranger()}.
#'
#' @return A fitted model object from \code{lm()} or \code{ranger()}.
#' @export
fit_model <- function(model, formula, synthetic_data, original_data, sensitive_attribute,
                      lm.control = list(), ranger.control = list()) {
  
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
                      probability = use_prob#,
                  #    na.action = na.omit
                    ),
                    ranger.control
                  )
                  do.call(ranger::ranger, args)
                },
                
                stop("Unsupported model type. Use 'lm' or 'rf'.")
  )
  
  return(fit)
}