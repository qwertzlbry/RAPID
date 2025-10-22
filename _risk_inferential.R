#' @title Inferential Disclosure Risk Assessment
#' 
#' @description
#' Assess inferential disclosure risk for a sensitive variable in a synthetic dataset. 
#' This function fits a model to predict a sensitive variable using known variables 
#' and measures how well the synthetic data allows an adversary to infer the true value. 
#' Supports continuous and categorical sensitive variables.
#'
#' @param truth The original dataset (data frame).
#' @param synth The synthetic dataset (data frame).
#' @param known_vars Character vector of quasi-identifiers used for inference.
#' @param sensitive_var Name of the sensitive variable to be predicted.
#' @param model Model type to use for inference: "lm" (linear) or "rf" (random forest).
#' @param method Error metric for continuous variables: "mae" or "rmse".
#' @param numeric_threshold_type Either "Percentage" or "Value".
#' @param numeric_threshold Threshold value (numeric).
#' @param categorical_threshold Threshold for categorical risk (numeric).
#' @param seed Optional random seed.
#' @param trace Logical; whether to print diagnostic messages.
#' @param ... Additional arguments passed to model fitting.
#' @author Oscar Thees, Matthias Templ
#'
#' @return A list containing:
#' - For continuous: number and percentage of low-error predictions and the subset of data with those predictions
#' - For categorical: confidence metrics for true class predictions
#'
#' @importFrom stats predict lm as.formula
#' @importFrom utils head
#' @importFrom ranger ranger
#' @examples
#' # Load required libraries
#' if (!requireNamespace("ranger", quietly = TRUE)) install.packages("ranger")
#' library(ranger)
#' 
#' # Generate toy original dataset (truth)
#' set.seed(2025)
#' n <- 100
#' truth <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, mean = 50000, sd = 10000)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   health_score = round(rnorm(n, mean = 75, sd = 10), 1)
#' )
#' 
#' # Generate a synthetic dataset with some variability
#' synth <- data.frame(
#'   age = truth$age + sample(-2:2, n, replace = TRUE),
#'   income = round(truth$income * runif(n, 0.9, 1.1)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   health_score = round(truth$health_score + rnorm(n, 0, 5), 1)
#' )
#' 
#' # Make sure factor levels are consistent
#' synth$education <- factor(synth$education, levels = levels(truth$education))
#' 
#' # Run the inferential disclosure risk function
#' result <- risk_inferential(
#'   truth = truth,
#'   synth = synth,
#'   known_vars = c("age", "income", "education"),
#'   sensitive_var = "health_score",
#'   model = "rf",
#'   method = "rmse",                      # root mean squared error
#'   numeric_threshold_type = "Percentage",
#'   numeric_threshold = 10,             # 10% relative error threshold
#'   trace = TRUE
#' ) # TODO: check the warning
#' 
#' # Output
#' print(result$rows_risk_n)
#' print(result$rows_risk_p)
#' head(result$rows_risk_df)
#' 
#' ## Categorical sensitive variable:
#' # Create toy truth data (original)
#' set.seed(2025)
#' n <- 100
#' truth_cat <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, 50000, 10000)),
#'   education = sample(c("low", "medium", "high"), n, replace = TRUE),
#'   disease_status = factor(sample(c("healthy", "diabetic", "hypertensive"), 
#'   n, replace = TRUE, prob = c(0.6, 0.2, 0.2)))
#' )
#' 
#' # Create toy synth data (synthetic)
#' synth_cat <- truth_cat
#' synth_cat$age <- synth_cat$age + sample(-2:2, n, replace = TRUE)
#' synth_cat$income <- round(synth_cat$income * runif(n, 0.9, 1.1))
#' synth_cat$education <- sample(c("low", "medium", "high"), n, replace = TRUE)
#' synth_cat$disease_status <- factor(sample(c("healthy", "diabetic", "hypertensive"), 
#' n, replace = TRUE, prob = c(0.55, 0.25, 0.2)))
#' 
#' # Apply the inferential_risk function with categorical sensitive variable
#' result_cat <- risk_inferential(
#'   truth = truth_cat,
#'   synth = synth_cat,
#'   known_vars = c("age", "income", "education"),
#'   sensitive_var = "disease_status",
#'   model = "rf",
#'   method = "rmse",
#'   numeric_threshold_type = "Percentage",
#'   numeric_threshold = 5,
#'   categorical_threshold = 1.2
#' )
#' 
#' # Print results
#' str(result_cat)
#' 
#' ## Compare old function (to be deleted)
#' result_cat <- inferential_risk(
#'   truth = truth_cat,
#'   synth = synth_cat,
#'   known_vars = c("age", "income", "education"),
#'   sensitive_var = "disease_status",
#'   model = "rf",
#'   method = "Euclidean",
#'   numeric_threshold_type = "Percentage",
#'   numeric_threshold = 5,
#'   categorical_threshold = 1.2
#' )
#' 
#' # Print results
#' str(result_cat)
#' 
#' @export
#' 
# Load all subfunctions, not needed anymore in a package environment. 
source("code/auxiliary_functions/risk/fit_model.R")
source("code/auxiliary_functions/risk/predict_sensitive.R")
source("code/auxiliary_functions/risk/compute_model_metrics.R")
source("code/auxiliary_functions/risk/evaluate_categorical_risk.R")
source("code/auxiliary_functions/risk/evaluate_continous_risk.R")
source("code/auxiliary_functions/risk/handle_sensitive_na.R")
source("code/auxiliary_functions/risk/handle_known_na.R")
source("code/auxiliary_functions/risk/ensure_all_levels.R")

risk_inferential <- function(truth,
                             synth,
                             known_vars,
                             sensitive_var,
                             numeric_na_strategy = c("constant", "drop", "median"),
                             constant_value = 0,
                             model = c("lm", "rf"),
                             method = c('mae', 'rmse', 'rmae', 'rrmse'),
                             numeric_threshold_type = c("Percentage", "Value"),
                             numeric_threshold = NULL,
                             categorical_threshold = 1,
                             seed = 2025,
                             trace = FALSE,
                             ...) {
  
  # Timing
  if(trace) {
    time_start <- Sys.time()
  }

  # Helper for message output
  log_msg <- function(...) if (trace) message(...)  
  
  # Argument matching and setup
  model <- match.arg(model)
  method <- tolower(method)
  method <- match.arg(method, choices = c("mae", "rmse", "rmae", "rrmse"))
  numeric_threshold_type <- match.arg(numeric_threshold_type)
  if (!is.null(seed)) set.seed(seed)
 
  stopifnot(is.data.frame(truth), is.data.frame(synth))
  stopifnot(all(c(known_vars, sensitive_var) %in% names(truth)))
  stopifnot(all(c(known_vars, sensitive_var) %in% names(synth)))
  
  # Compatibility check: lm is not appropriate for categorical targets
  if (model == "lm" && is.factor(truth[[sensitive_var]])) {
    stop("Model 'lm' is not suitable for categorical sensitive variables. Use 'rf' or another classification model instead.")
  }
  
  # Subset to needed variables--------------------------------------------------
  truth <- truth[, c(sensitive_var, known_vars)]
  synth <- synth[, c(sensitive_var, known_vars)]
  
  # Handle NA in sensitive and known variables (sub function)-------------------
  res <- handle_sensitive_na(truth, synth, sensitive_var, numeric_na_strategy, constant_value)
  truth <- res$truth
  synth <- res$synth
  
  # only  replacement with NA catergoies
  synth <- handle_known_na(synth, known_vars)
  truth <- handle_known_na(truth, known_vars)
  
  # ensure all levels are used for fitting
  synth <- ensure_all_levels_used(synth, truth, sensitive_var)

  # Fit model (sub function)----------------------------------------------------
  
  ## Formula
  formula <- as.formula(paste(sensitive_var, "~", paste(known_vars, collapse = "+")))
  
  fit <- fit_model(model, formula, synth, truth, sensitive_var,
                        lm.control = list(), ranger.control = list()) 
  
  # Predict (sub function)------------------------------------------------------
  A <- truth[[sensitive_var]]
  B <- predict_sensitive(model, fit, truth, sensitive_var)
  
  # print("Check predicted class matrix B:")
  # print(colnames(B))
  # print("Unique A:")
  # print(unique(as.character(A)))
  # print("Rows with A == 'missing':")
  # print(sum(as.character(A) == "missing"))
  # print("Mean of B[, 'missing'] for those rows:")
  # print(mean(B[as.character(A) == "missing", "missing"]))
  
  # Model metrics (sub function) #provisional-----------------------------------
  metrics <- compute_model_metrics(A, B)
  
  # Categorical evaluation (sub function)---------------------------------------
  if (is.factor(A)) {
    eval <- evaluate_categorical(A, B, categorical_threshold)
    #  comment: for categorical sensitive variables, put additional output
    #           (see function interential_risk()
  } else {
    
  # Continuous evaluation (sub function)----------------------------------------
    eval <- evaluate_numeric(A, B, truth, method, numeric_threshold, numeric_threshold_type)
  }
  
  if (trace){
    message("Time: ", round(difftime(Sys.time(), time_start, 
                                     units = "secs"), 2), " seconds")
  }
  return(list(
    risk = eval,         # output from evaluate_numeric() or evaluate_categorical_prediction()
    metrics = metrics     # output from compute_model_metrics()
  ))
} # End of function
