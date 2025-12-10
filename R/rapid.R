#' @title RAPID: Risk of Attribute Prediction–Induced Disclosure
#' @author Oscar Thees, Matthias Templ
#'
#' @description
#' Assess inferential disclosure risk for a sensitive variable in a synthetic dataset using the RAPID metric.
#' This function trains a predictive model on synthetic data using quasi-identifiers, then evaluates attribute inference risk
#' by scoring the model on the original data. Supports continuous and categorical sensitive variables.
#'
#' @param original_data The original dataset (data frame).
#' @param synthetic_data The synthetic dataset (data frame).
#' @param quasi_identifiers Character vector of quasi-identifiers used for inference.
#' @param sensitive_attribute Name of the sensitive variable to be predicted.
#' @param model_type Model type to use for inference: "lm" (linear) or "rf" (random forest).
#' @param num_error_metric Error metric for continuous variables: "mae", "rmse", "rmae", or "rrmse".
#' @param num_epsilon_type Either "Percentage" or "Value".
#' @param num_epsilon Threshold for continuous attributes (numeric).
#' @param cat_tau Threshold for categorical risk (numeric).
#' @param seed Optional random seed.
#' @param trace Logical; whether to print diagnostic messages.
#' @param ... Additional arguments passed to model fitting.
#'
#'
#' @return A list containing:
#' - For continuous: proportion and subset of records within error tolerance.
#' - For categorical: RAPID metric and supporting confidence measures.
#'
#'
#' @importFrom stats predict lm as.formula
#' @importFrom utils head
#' @importFrom ranger ranger
#' @examples
#' # Load required libraries
#' if (!requireNamespace("ranger", quietly = TRUE)) install.packages("ranger")
#' library(ranger)
#'
#' # Generate toy original dataset
#' set.seed(2025)
#' n <- 100
#' original_data <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, mean = 50000, sd = 10000)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   health_score = round(rnorm(n, mean = 75, sd = 10), 1)
#' )
#'
#' # Generate synthetic dataset with noise
#' synthetic_data <- data.frame(
#'   age = original_data$age + sample(-2:2, n, replace = TRUE),
#'   income = round(original_data$income * runif(n, 0.9, 1.1)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   health_score = round(original_data$health_score + rnorm(n, 0, 5), 1)
#' )
#'
#' # Match factor levels
#' synthetic_data$education <- factor(synthetic_data$education, levels = levels(original_data$education))
#'
#' # Run RAPID for continuous sensitive attribute
#' result <- rapid(
#'   original_data = original_data,
#'   synthetic_data = synthetic_data,
#'   quasi_identifiers = c("age", "income", "education"),
#'   sensitive_attribute = "health_score",
#'   model_type = "rf",
#'   num_error_metric = "rmse",                  # root mean squared error
#'   num_epsilon_type = "Percentage",
#'   num_epsilon = 10,                           # 10% relative error threshold
#'   trace = TRUE
#' )
#'
#' # Output
#' print(result$risk$rows_risk_n)
#' print(result$risk$rows_risk_p)
#' head(result$risk$rows_risk_df)
#'
#' # ---- Categorical sensitive variable example ----
#' set.seed(2025)
#' original_data_cat <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, 50000, 10000)),
#'   education = sample(c("low", "medium", "high"), n, replace = TRUE),
#'   disease_status = factor(sample(c("healthy", "diabetic", "hypertensive"),
#'     n, replace = TRUE, prob = c(0.6, 0.2, 0.2)))
#' )
#'
#' synthetic_data_cat <- original_data_cat
#' synthetic_data_cat$age <- synthetic_data_cat$age + sample(-2:2, n, replace = TRUE)
#' synthetic_data_cat$income <- round(synthetic_data_cat$income * runif(n, 0.9, 1.1))
#' synthetic_data_cat$education <- sample(c("low", "medium", "high"), n, replace = TRUE)
#' synthetic_data_cat$disease_status <- factor(sample(c("healthy", "diabetic", "hypertensive"),
#'     n, replace = TRUE, prob = c(0.55, 0.25, 0.2)))
#'
#' # Run RAPID for categorical sensitive attribute
#' result_cat <- rapid(
#'   original_data = original_data_cat,
#'   synthetic_data = synthetic_data_cat,
#'   quasi_identifiers = c("age", "income", "education"),
#'   sensitive_attribute = "disease_status",
#'   model_type = "rf",
#'   cat_tau = 1.2
#' )
#'
#' # Print result summary
#' str(result_cat)
#'
#' @export
rapid <- function(original_data,
                  synthetic_data,
                  quasi_identifiers,
                  sensitive_attribute,
                  model_type = c("lm", "rf", "cart", "gbm", "logit"),

                  # Numeric-specific
                  num_na_strategy = c("constant", "drop", "median"), #num_na_strategy
                  num_constant_value = 0, #num_constant_value
                  num_error_metric = c("mae", "rmse", "rmae", "rrmse"), #num_error_metric
                  num_epsilon_type = c("Percentage", "Value"), # num_epsilon_type
                  num_epsilon = NULL, #num_epsilon

                  # Categorical-specific
                  cat_eval_method =c("RCS_conditional", "RCS_marginal", "NCE"), #cat_eval_method
                  cat_tau = 1, #cat_tau

                  # Varia
                  seed = 2025,
                  trace = FALSE,
                  return_all_records = FALSE,
                  ...) {

  if (trace){
  time_start <- Sys.time()
  }

  log_msg <- function(...) if (trace) message(...)

  # Argument matching and setup
  model_type <- match.arg(model_type)
  num_error_metric <- tolower(num_error_metric)
  num_error_metric <- match.arg(num_error_metric, choices = c("mae", "rmse", "rmae", "rrmse"))
  num_epsilon_type <- match.arg(num_epsilon_type)
  num_na_strategy <- match.arg(num_na_strategy)
  cat_eval_method <- match.arg(cat_eval_method)


  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Validation
  stopifnot(is.data.frame(original_data), is.data.frame(synthetic_data))
  stopifnot(all(c(quasi_identifiers, sensitive_attribute) %in% names(original_data)))
  stopifnot(all(c(quasi_identifiers, sensitive_attribute) %in% names(synthetic_data)))

  sensitive_attr_vec <- original_data[[sensitive_attribute]]
  is_categorical <- is.factor(sensitive_attr_vec) || is.character(sensitive_attr_vec)


  # Compatibility check
  if (model_type == "lm" && is.factor(original_data[[sensitive_attribute]])) {
    stop("Model 'lm' is not suitable for categorical sensitive variables. Use 'rf'
         or another classification model instead.")
  }
  if (trace) {
      cat("  RAPID: Running Single Evaluation Mode\n")
    }
  # Subset to needed variables--------------------------------------------------
  original_data <- original_data[, c(sensitive_attribute, quasi_identifiers)]
  synthetic_data <- synthetic_data[, c(sensitive_attribute, quasi_identifiers)]

  # Handle NA in sensitive and known variables (sub function)-------------------
  res <- handle_sensitive_na(original_data, synthetic_data, sensitive_attribute,
                             num_na_strategy, num_constant_value)
  original_data <- res$original_data
  synthetic_data <- res$synthetic_data

  # only  replacement with NA catergoies
  synthetic_data <- handle_known_na(synthetic_data, quasi_identifiers)
  original_data <- handle_known_na(original_data, quasi_identifiers)

  # ensure all levels are used for fitting
  synthetic_data <- ensure_all_levels_used(synthetic_data, original_data, sensitive_attribute)

  # Fit model (sub function)----------------------------------------------------
  formula <- as.formula(paste(sensitive_attribute, "~", paste(quasi_identifiers, collapse = "+")))
  fit <- fit_model(model_type,
                   formula, synthetic_data,
                   original_data,
                   sensitive_attribute,
                   lm.control = list(),
                   ranger.control = list(),
                   rpart.control = list(),
                   xgb.control = list()
                   )

   # Predict on original data (sub function)------------------------------------
  A <- original_data[[sensitive_attribute]]
  B <- predict_sensitive(model_type, fit, original_data, sensitive_attribute)

  # Model metrics (sub function) #provisional-----------------------------------
  metrics <- compute_model_metrics(A, B)

  # Categorical evaluation (sub function)---------------------------------------
  if (is.factor(A)) {
    eval <- evaluate_categorical(A, B, cat_tau, original_data, cat_eval_method, sensitive_attribute, return_all_records)
  } else {
  # Continuous evaluation (sub function)----------------------------------------
    eval <- evaluate_numeric(A, B, original_data, num_error_metric, num_epsilon, num_epsilon_type, return_all_records)
  }

  if (trace){
    message("Time: ", round(difftime(Sys.time(), time_start, units = "secs"), 2), " seconds")
  }
  return(list(
    risk = eval, # output from evaluate_numeric() or evaluate_categorical_prediction()
    metrics = metrics  # output from compute_model_metrics()
  ))
}

#' @rdname rapid
#' @export
risk_inferential <- function(...) rapid(...)

