#' @title RAPID: Risk of Attribute Prediction-Induced Disclosure
#' @author Oscar Thees, Matthias Templ
#'
#' @description
#' Assess inferential disclosure risk for a sensitive variable in a synthetic dataset using the RAPID metric.
#' The function trains an attacker model on synthetic data using quasi-identifiers and evaluates attribute inference risk
#' by scoring the attacker on the original data. Supports continuous and categorical sensitive variables.
#'
#' @param original_data The original dataset (data.frame).
#' @param synthetic_data The synthetic dataset (data.frame).
#' @param quasi_identifiers Character vector of quasi-identifiers used for inference.
#' @param sensitive_attribute Name (character scalar) of the sensitive variable to be predicted.
#' @param model_type Model type for the attacker:
#'   \code{"rf"} (random forest, default),
#'   \code{"cart"} (classification/regression tree),
#'   \code{"gbm"} (gradient boosting machine),
#'   \code{"logit"} (logistic regression; categorical only),
#'   \code{"lm"} (linear model; continuous only).
#'
#' @param num_na_strategy Strategy for handling missing values in a continuous sensitive attribute:
#'   \code{"constant"} (replace with a constant),
#'   \code{"drop"} (remove affected records),
#'   \code{"median"} (impute with the median).
#' @param num_constant_value Constant used when \code{num_na_strategy = "constant"}. Default is 0.
#'
#' @param num_error_metric Error metric for continuous sensitive attributes:
#'   \code{"symmetric"} (recommended; symmetric relative error),
#'   \code{"stabilised_relative"} (stabilised relative error),
#'   \code{"absolute"} (absolute error on the original scale).
#' @param num_epsilon Numeric threshold for continuous attributes. Interpretation depends on \code{num_error_metric}:
#'   \itemize{
#'     \item For \code{"symmetric"} and \code{"stabilised_relative"}, \code{num_epsilon} is a relative error threshold
#'     on a 0--1 scale (e.g., \code{0.05} corresponds to 5\%).
#'     \item For \code{"absolute"}, \code{num_epsilon} is an absolute error threshold on the original scale
#'     of the sensitive attribute.
#'   }
#' @param num_delta Smoothing constant to prevent division by zero in relative error metrics. Default is 0.01.
#'
#' @param cat_eval_method Evaluation method for categorical sensitive attributes:
#'   \code{"RCS_marginal"} (recommended; normalized gain over marginal baseline),
#'   \code{"RCS_conditional"} (ratio over class-conditional baseline),
#'   \code{"NCE"} (normalized cross-entropy).
#' @param cat_tau Numeric threshold for categorical risk. Interpretation depends on \code{cat_eval_method}.
#'
#' @param seed Optional random seed for reproducibility.
#' @param trace Logical; if \code{TRUE}, print diagnostic messages during execution.
#' @param return_all_records Logical; if \code{TRUE}, return all records with their risk status in \code{rows_risk_df}.
#'   If \code{FALSE} (default), return only at-risk records.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with two components:
#' \describe{
#'   \item{risk}{Output from \code{evaluate_numeric()} or \code{evaluate_categorical()}, containing (at least):
#'     \itemize{
#'       \item \code{confidence_rate}: Proportion of records at risk (0--1)
#'       \item \code{n_at_risk}: Number of records flagged as at-risk
#'       \item \code{percentage}: Percentage of records at risk (0--100)
#'       \item \code{rows_risk_df}: Data frame of records with risk assessments (all or at-risk only)
#'     }
#'   }
#'   \item{metrics}{Model performance metrics from \code{compute_model_metrics()}, e.g. accuracy/F1 for categorical
#'     or RMSE/MAE for continuous outcomes (depending on implementation).}
#' }
#'
#' @details
#' RAPID (Risk of Attribute Prediction-Induced Disclosure) quantifies the risk that an attacker can infer sensitive
#' attributes from released synthetic data by training a predictive model on the synthetic data.
#' The attacker is evaluated on the original data; records are flagged as at-risk when the prediction error
#' is below a user-chosen threshold.
#'
#' Categorical sensitive attributes can be evaluated via \code{RCS_marginal}, \code{RCS_conditional}, or \code{NCE}.
#' Continuous sensitive attributes can be evaluated via relative error metrics (\code{"symmetric"},
#' \code{"stabilised_relative"}) or absolute error (\code{"absolute"}).
#'
#' @importFrom stats predict lm as.formula
#' @importFrom utils head
#' @importFrom ranger ranger
#'
#' @examples
#' \dontrun{
#' library(ranger)
#' library(synthpop)
#'
#' set.seed(2025)
#' n <- 100
#' original_data <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, mean = 50000, sd = 10000)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   health_score = round(rnorm(n, mean = 75, sd = 10), 1)
#' )
#'
#' synthetic_data <- syn(original_data, method = "cart")$syn
#'
#' # ---- Continuous sensitive attribute (relative error) ----
#' result_cont <- rapid(
#'   original_data = original_data,
#'   synthetic_data = synthetic_data,
#'   quasi_identifiers = c("age", "income", "education"),
#'   sensitive_attribute = "health_score",
#'   model_type = "rf",
#'   num_error_metric = "symmetric",
#'   num_epsilon = 0.05,  # 5% threshold (0--1 scale)
#'   return_all_records = TRUE
#' )
#'
#' print(result_cont$risk$confidence_rate)
#' head(result_cont$risk$rows_risk_df)
#'
#' # ---- Continuous sensitive attribute (absolute error) ----
#' result_abs <- rapid(
#'   original_data = original_data,
#'   synthetic_data = synthetic_data,
#'   quasi_identifiers = c("age", "income", "education"),
#'   sensitive_attribute = "health_score",
#'   model_type = "rf",
#'   num_error_metric = "absolute",
#'   num_epsilon = 2,      # absolute units of health_score
#'   return_all_records = TRUE
#' )
#'
#' # ---- Categorical sensitive attribute ----
#' original_data_cat <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = round(rnorm(n, 50000, 10000)),
#'   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#'   disease_status = factor(sample(c("healthy", "diabetic", "hypertensive"),
#'     n, replace = TRUE, prob = c(0.6, 0.2, 0.2)))
#' )
#'
#' synthetic_data_cat <- syn(original_data_cat, method = "cart")$syn
#'
#' result_cat <- rapid(
#'   original_data = original_data_cat,
#'   synthetic_data = synthetic_data_cat,
#'   quasi_identifiers = c("age", "income", "education"),
#'   sensitive_attribute = "disease_status",
#'   model_type = "rf",
#'   cat_eval_method = "RCS_marginal",
#'   cat_tau = 0.3,
#'   return_all_records = TRUE
#' )
#'
#' print(result_cat$risk$confidence_rate)
#' head(result_cat$risk$rows_risk_df)
#' }
#'
#' @export

rapid <- function(original_data,
                  synthetic_data,
                  quasi_identifiers,
                  sensitive_attribute,
                  model_type = c("rf", "cart", "gbm", "logit","lm"),

                  # Numeric-specific
                  num_na_strategy = c("constant", "drop", "median"), #num_na_strategy
                  num_constant_value = 0, #num_constant_value
                  num_error_metric = c("symmetric", "stabilised_relative", "absolute"),
                  num_epsilon = NULL, #num_epsilon
                  num_epsilon_scale = c("proportion", "percent"),
                  num_delta = 0.01,

                  # Categorical-specific
                  cat_eval_method =c("RCS_marginal","RCS_conditional", "NCE"), #cat_eval_method
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
  num_error_metric <- match.arg(num_error_metric)
  num_epsilon_scale <- match.arg(num_epsilon_scale)
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
    eval <- evaluate_numeric(A, B, original_data, num_error_metric, num_epsilon, num_epsilon_scale, num_delta, return_all_records)
  }

  if (trace){
    message("Time: ", round(difftime(Sys.time(), time_start, units = "secs"), 2), " seconds")
  }

  result <- list(
    risk = eval, # output from evaluate_numeric() or evaluate_categorical_prediction()
    metrics = metrics  # output from compute_model_metrics()
  )
  class(result) <- "rapid_result"
return(result)
}

#' @rdname rapid
#' @export
risk_inferential <- function(...) rapid(...)

