
################################################################################

# library(riskutility)   # Your package with RAPID implementation
# Load required packages

devtools::load_all()
library(synthpop)      # For synthetic data generation (CART)
library(ranger)        # For Random Forest
library(boot)          # For bootstrapped confidence intervals
library(dplyr)         # For data wrangling
library(ggplot2)       # For plotting
library(caret)         # For calibration plots

# Load UCI Adult dataset
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
df <- read.csv(url, header = FALSE, strip.white = TRUE)

# Set column names as per UCI documentation
colnames(df) <- c(
  "age", "workclass", "fnlwgt", "education", "education.num", "marital.status",
  "occupation", "relationship", "race", "sex", "capital.gain", "capital.loss",
  "hours.per.week", "native.country", "income"
)

# Convert character columns to factors (nominal scale)
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], factor)

# Drop unneeded or redundant variables
df <- df %>% select(-fnlwgt, -education)

# Remove rows with missing values in the sensitive attribute
df <- df[complete.cases(df$income), ]

# Generate synthetic datasets using CART
synths <- synthpop::syn(df, m = 5, method = "cart")
synthetic_list <- synths$syn  # List of 5 synthetic datasets

# Define function to compute RAPID for one replicate
compute_rapid_cat <- function(original_data, synthetic_data, y = "income", cat_tau = 0.3, model = "rf", n_boot = 1000) {
  # Fit model on synthetic data and evaluate on original
  rapid_result <- rapid(
    original_data = original_data,
    synthetic_data = synthetic_data,
    quasi_identifiers = setdiff(colnames(original_data), y),
    sensitive_attribute = y,
    model_type = model,
    cat_tau = cat_tau,
    trace = FALSE
  )

  # Bootstrap RCIR
  rapid_vals <- boot::boot(
    data = original_data,
    statistic = function(data, i) {
      tmp_orig <- data[i, ]
      rapid_i <- rapid(
        original_data = tmp_orig,
        synthetic_data = synthetic_data,
        quasi_identifiers = setdiff(colnames(original_data), y),
        sensitive_attribute = y,
        model_type = model,
        cat_tau = cat_tau,
        trace = FALSE
      )
      return(rapid_i$risk$confidence_rate)
    },
    R = n_boot
  )

  ci <- boot.ci(rapid_vals, type = "perc")$percent[4:5]

  return(list(
    rapid = rapid_result$risk$confidence_rate,
    ci = ci
  ))
}

# Run across 5 synthetic replicates
rapid_results <- lapply(synthetic_list, function(syn_df) {
  compute_rapid_cat(original_data = df, synthetic_data = syn_df, cat_tau = 0.3, n_boot = 1000)
})

# Aggregate and print results
rapid_summary <- data.frame(
  replicate = 1:5,
  RCIR = sapply(rapid_results, function(x) x$rapid),
  CI_Low = sapply(rapid_results, function(x) x$ci[1]),
  CI_High = sapply(rapid_results, function(x) x$ci[2])
)

print(rapid_summary)

# Plot RCIR threshold curve (τ)
taus <- seq(0.05, 1.0, by = 0.05)

# Compute RAPID for each tau
curve_vals <- sapply(taus, function(t) {
  compute_rapid_cat(
    original_data = df,
    synthetic_data = synthetic_list[[1]],
    y = "income",
    cat_tau = t,
    n_boot = 100
  )$rapid
})

pdf("/Users/matthias/Downloads/rcir_threshold_curve.pdf", width = 8, height = 6)
ggplot(data.frame(tau = taus, rapid = curve_vals), aes(x = tau, y = rapid)) +
  geom_line() + geom_point() +
  #geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", size = 1)
  labs(title = "RCIR vs Confidence Threshold (τ)",
       x = "Confidence threshold τ", y = "RCIR") +
  theme_minimal()
dev.off()


# Compute RCIR for each tau and each synthetic replicate
rapid_matrix <- sapply(taus, function(tau_val) {
  sapply(synthetic_list, function(synth_df) {
    compute_rapid_cat(
      original_data = df,
      synthetic_data = synth_df,
      y = "income",
      tau = tau_val,
      n_boot = 10  # keep small for demo; increase to 500 for final
    )$rapid
  })
})

# Transpose so rows = tau, columns = replicates
rapid_df <- as.data.frame(t(rcir_matrix))
colnames(rcir_df) <- paste0("rep", 1:ncol(rapid_df))
rapid_df$tau <- taus

# Compute mean, 5% and 95% quantiles across replicates
library(dplyr)
rapid_summary <- rapid_df %>%
  tidyr::pivot_longer(cols = starts_with("rep"), names_to = "replicate", values_to = "rcir") %>%
  group_by(tau) %>%
  summarise(
    rapid_mean = mean(rapid),
    rapid_q05  = quantile(rapid, 0.05),
    rapid_q95  = quantile(rapid, 0.95),
    .groups = "drop"
  )

# Plot mean curve with quantile ribbon
library(ggplot2)
pdf("/Users/matthias/Downloads/rcir_curve_quantile.pdf", width = 8, height = 6)
ggplot(rapid_summary, aes(x = tau, y = rapid_mean)) +
  geom_line(color = "black", size = 0.3) +
  geom_ribbon(aes(ymin = rapid_q05, ymax = rapid_q95),
              fill = "darkorange", alpha = 0.3) +
  geom_point(size = 1.2) +
  labs(
    title = expression(RAPID ~ "vs Confidence Threshold (" * tau * ")"),
    subtitle = "Mean ± Confidence Interval (Quantile Range (0.05–0.95)) across 5 synthetic replicates",
    x = expression(paste("Confidence Threshold ", tau)),
    y = "RCIR"
  ) +
  theme_minimal(base_size = 13)
dev.off()



################################################################################
###   Table 1   ####
###############################################################################

rm(list = ls())
sapply(list.files("code/auxiliary_functions/risk_paper", pattern = "\\.R$", full.names = TRUE), source)
library(dplyr)
library(synthpop)
library(boot)
library(riskutility)

# Load UCI Adult dataset
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
df <- read.csv(url, header = FALSE, strip.white = TRUE)

# Set column names as per UCI documentation
colnames(df) <- c(
  "age", "workclass", "fnlwgt", "education", "education.num", "marital.status",
  "occupation", "relationship", "race", "sex", "capital.gain", "capital.loss",
  "hours.per.week", "native.country", "income"
)

# Convert character columns to factors (nominal scale)
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], factor)

# Drop unneeded or redundant variables
df <- df %>% select(-fnlwgt, -education)

# Remove rows with missing values in the sensitive attribute
df <- df[complete.cases(df$income), ]

# Generate synthetic datasets using CART
synths <- synthpop::syn(df, m = 5, method = "cart")
synthetic_list <- synths$syn  # List of 5 synthetic datasets


# Assuming: compute_rapid_cat() is already defined

# Define attacker models and thresholds
attackers <- c("cart","rf", "gbm", "logit")  # Names used in your rapid() function
tau_vals <- c(1.10, 1.25, 1.50)
n_boot <- 10  # Bootstrap replicates

# Create results container
results <- list()


rapid <- function(original_data,
                  synthetic_data,
                  quasi_identifiers,
                  sensitive_attribute,
                  numeric_na_strategy = c("constant", "drop", "median"),
                  constant_value = 0,
                  model_type = c("lm", "rf", "gbm", "logit", "cart"),
                  error_metric = c("mae", "rmse", "rmae", "rrmse"),
                  epsilon_type = c("Percentage", "Value"),
                  epsilon = NULL,
                  tau = 1,
                  seed = 2025,
                  trace = FALSE,
                  ...) {

  if(trace) {
    time_start <- Sys.time()
  }

  log_msg <- function(...) if (trace) message(...)

  model_type <- match.arg(model_type)
  error_metric <- tolower(error_metric)
  error_metric <- match.arg(error_metric, choices = c("mae", "rmse", "rmae", "rrmse"))
  epsilon_type <- match.arg(epsilon_type)
  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.data.frame(original_data), is.data.frame(synthetic_data))
  stopifnot(all(c(quasi_identifiers, sensitive_attribute) %in% names(original_data)))
  stopifnot(all(c(quasi_identifiers, sensitive_attribute) %in% names(synthetic_data)))

  # Disallow 'lm' for categorical response
  if (model_type == "lm" && is.factor(original_data[[sensitive_attribute]])) {
    stop("Model 'lm' is not suitable for categorical sensitive variables. Use 'rf', 'gbm', 'logit', or 'cart'.")
  }

  # Prepare data
  original_data <- original_data[, c(sensitive_attribute, quasi_identifiers)]
  synthetic_data <- synthetic_data[, c(sensitive_attribute, quasi_identifiers)]

  res <- handle_sensitive_na(original_data, synthetic_data, sensitive_attribute, numeric_na_strategy, constant_value)
  original_data <- res$original_data
  synthetic_data <- res$synthetic_data

  synthetic_data <- handle_known_na(synthetic_data, quasi_identifiers)
  original_data <- handle_known_na(original_data, quasi_identifiers)
  synthetic_data <- ensure_all_levels_used(synthetic_data, original_data, sensitive_attribute)

  res2 <- align_factor_levels_for_glm(original_data, synthetic_data, quasi_identifiers)
  original_data <- res2$original_data
  synthetic_data <- res2$synthetic_data

  formula <- as.formula(paste(sensitive_attribute, "~", paste(quasi_identifiers, collapse = "+")))

  # Fit attacker model
  fit <- fit_model(model_type, formula, synthetic_data, original_data, sensitive_attribute,
                   lm.control = list(), ranger.control = list(), ...)

  A <- original_data[[sensitive_attribute]]
  B <- predict_sensitive(model_type, fit, original_data, sensitive_attribute)

  metrics <- compute_model_metrics(A, B)

  if (is.factor(A)) {
    eval <- evaluate_categorical(A, B, tau)
  } else {
    eval <- evaluate_numeric(A, B, original_data, error_metric, epsilon, epsilon_type)
  }

  if (trace){
    message("Time: ", round(difftime(Sys.time(), time_start, units = "secs"), 2), " seconds")
  }

  return(list(
    risk = eval,
    metrics = metrics
  ))
}

inject_missing_factor_levels <- function(train_data, reference_data, vars) {
  for (v in vars) {
    if (is.factor(train_data[[v]]) && is.factor(reference_data[[v]])) {
      all_levels <- union(levels(train_data[[v]]), levels(reference_data[[v]]))
      levels(train_data[[v]]) <- all_levels

      # Add one row for each level not present in training data
      missing_levels <- setdiff(levels(reference_data[[v]]), unique(train_data[[v]]))
      for (lvl in missing_levels) {
        ghost_row <- train_data[1, , drop = FALSE]
        ghost_row[[v]] <- factor(lvl, levels = all_levels)
        train_data <- rbind(train_data, ghost_row)
      }
    }
  }
  return(train_data)
}

align_factor_levels_for_glm <- function(original_data, synthetic_data, vars) {
  for (v in vars) {
    if (is.factor(original_data[[v]]) && is.factor(synthetic_data[[v]])) {
      all_levels <- union(levels(original_data[[v]]), levels(synthetic_data[[v]]))

      # Add missing levels to synthetic_data *before model fitting*
      levels(synthetic_data[[v]]) <- all_levels

      # Optionally harmonize original_data as well
      levels(original_data[[v]]) <- all_levels
    }
  }
  list(original_data = original_data, synthetic_data = synthetic_data)
}

# fit_model <- function(model_type, formula, train_data, test_data, response,
#                       lm.control = list(), ranger.control = list(), ...) {
#
#   if (model_type == "lm") {
#     return(do.call(stats::lm, c(list(formula = formula, data = train_data), lm.control)))
#
#   } else if (model_type == "rf") {
#     library(ranger)
#     return(do.call(ranger::ranger, c(list(formula = formula, data = train_data, probability = TRUE), ranger.control)))
#
#   } else if (model_type == "gbm") {
#     library(xgboost)
#
#     # Prepare data for xgboost
#     X <- model.matrix(formula, data = train_data)[, -1]  # remove intercept
#     y <- as.numeric(train_data[[response]]) - 1           # convert factor to 0/1
#
#     dtrain <- xgb.DMatrix(data = X, label = y)
#
#     params <- list(
#       objective = "binary:logistic",
#       eval_metric = "logloss",
#       verbosity = 0
#     )
#
#     model <- xgboost::xgb.train(
#       params = params,
#       data = dtrain,
#       nrounds = 100,
#       verbose = 0
#     )
#
#     model$.__x_data__ <- X  # Attach for prediction function
#     model$.__x_formula__ <- formula
#     return(model)
#
#   } else if (model_type == "logit") {
#     # Inject ghost rows for missing levels
#     train_data_expanded <- inject_missing_factor_levels(train_data, test_data, quasi_identifiers)
#
#     # Fit model on expanded data
#     return(glm(formula = formula, data = train_data_expanded, family = binomial()))
#   } else {
#     stop("Unsupported model type. Use 'lm', 'rf', 'gbm', or 'logit'.")
#   }
# }

fit_model <- function(model_type, formula, train_data, test_data, response,
                      lm.control = list(), ranger.control = list(), ...) {

  if (model_type == "lm") {
    return(do.call(stats::lm, c(list(formula = formula, data = train_data), lm.control)))

  } else if (model_type == "rf") {
    library(ranger)
    return(do.call(ranger::ranger, c(list(formula = formula, data = train_data, probability = TRUE), ranger.control)))

  } else if (model_type == "gbm") {
    library(xgboost)

    X <- model.matrix(formula, data = train_data)[, -1]
    y <- as.numeric(train_data[[response]]) - 1

    dtrain <- xgb.DMatrix(data = X, label = y)

    params <- list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      verbosity = 0
    )

    model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = 100,
      verbose = 0
    )

    model$.__x_data__ <- X
    model$.__x_formula__ <- formula
    return(model)

  } else if (model_type == "logit") {
    train_data_expanded <- inject_missing_factor_levels(train_data, test_data, quasi_identifiers)
    return(glm(formula = formula, data = train_data_expanded, family = binomial()))

  } else if (model_type == "cart") {
    library(rpart)

    # Select method based on response type
    method_type <- if (is.factor(train_data[[response]])) "class" else "anova"

    control <- rpart::rpart.control(...)  # Allow passing control args via ...

    return(rpart::rpart(
      formula = formula,
      data = train_data,
      method = method_type,
      control = control
    ))

  } else {
    stop("Unsupported model type. Use 'lm', 'rf', 'gbm', 'logit', or 'cart'.")
  }
}

# predict_sensitive <- function(model_type, model, test_data, response) {
#   if (model_type == "lm") {
#     return(predict(model, newdata = test_data))
#
#   } else if (model_type == "rf") {
#     pred <- predict(model, data = test_data)$predictions
#     return(as.data.frame(pred))
#
#   } else if (model_type == "gbm") {
#     X_test <- model.matrix(model$.__x_formula__, data = test_data)[, -1]
#     probs <- predict(model, newdata = X_test)
#     probs_df <- data.frame(
#       `<=50K` = 1 - probs,
#       `>50K` = probs
#     )
#     return(probs_df)
#
#   } else if (model_type == "logit") {
#     probs <- predict(model, newdata = test_data, type = "response")
#     probs_df <- data.frame(
#       `<=50K` = 1 - probs,
#       `>50K` = probs
#     )
#     return(probs_df)
#
#   } else {
#     stop("Unsupported model type in predict_sensitive().")
#   }
# }

# predict_sensitive <- function(model_type, model, test_data, response) {
#   if (model_type == "lm") {
#     return(predict(model, newdata = test_data))
#
#   } else if (model_type == "rf") {
#     pred <- predict(model, data = test_data)$predictions
#     return(as.data.frame(pred))
#
#   } else if (model_type == "gbm") {
#     X_test <- model.matrix(model$.__x_formula__, data = test_data)[, -1]
#     probs <- predict(model, newdata = X_test)
#
#     # Get levels of original response
#     levels_response <- levels(test_data[[response]])
#     stopifnot(length(levels_response) == 2)
#
#     # Construct a valid probability matrix (2 columns: class0, class1)
#     probs_df <- data.frame(
#       prob0 = 1 - probs,
#       prob1 = probs
#     )
#     colnames(probs_df) <- levels_response  # Assign factor levels as column names
#
#     return(probs_df)
#
#   } else if (model_type == "logit") {
#     probs <- predict(model, newdata = test_data, type = "response")
#     levels_response <- levels(test_data[[response]])
#     stopifnot(length(levels_response) == 2)
#
#     probs_df <- data.frame(
#       prob0 = 1 - probs,
#       prob1 = probs
#     )
#     colnames(probs_df) <- levels_response
#     return(probs_df)
#
#   } else {
#     stop("Unsupported model type in predict_sensitive().")
#   }
# }

predict_sensitive <- function(model_type, model, test_data, response) {
  if (model_type == "lm") {
    return(predict(model, newdata = test_data))

  } else if (model_type == "rf") {
    pred <- predict(model, data = test_data)$predictions
    return(as.data.frame(pred))

  } else if (model_type == "gbm") {
    X_test <- model.matrix(model$.__x_formula__, data = test_data)[, -1]
    probs <- predict(model, newdata = X_test)

    levels_response <- levels(test_data[[response]])
    stopifnot(length(levels_response) == 2)

    probs_df <- data.frame(
      prob0 = 1 - probs,
      prob1 = probs
    )
    colnames(probs_df) <- levels_response
    return(probs_df)

  } else if (model_type == "logit") {
    probs <- predict(model, newdata = test_data, type = "response")
    levels_response <- levels(test_data[[response]])
    stopifnot(length(levels_response) == 2)

    probs_df <- data.frame(
      prob0 = 1 - probs,
      prob1 = probs
    )
    colnames(probs_df) <- levels_response
    return(probs_df)

  } else if (model_type == "cart") {
    # Use type = "prob" for rpart
    probs <- predict(model, newdata = test_data, type = "prob")

    # rpart returns a matrix with columns as levels
    probs_df <- as.data.frame(probs)

    # Ensure column names match response levels
    levels_response <- levels(test_data[[response]])
    if (!all(levels_response %in% colnames(probs_df))) {
      stop("CART model prediction does not contain all levels of the response.")
    }

    # Reorder columns to match factor levels
    probs_df <- probs_df[, levels_response, drop = FALSE]

    return(probs_df)

  } else {
    stop("Unsupported model type in predict_sensitive().")
  }
}

compute_model_metrics <- function(A, B) {
  # If B is a data frame of predicted probabilities, convert to predicted class
  if (is.data.frame(B)) {
    pred_class <- colnames(B)[apply(B, 1, which.max)]
    pred_class <- factor(pred_class, levels = levels(A))
  } else {
    pred_class <- B
  }

  acc <- mean(pred_class == A)
  return(list(accuracy = acc))
}

# Loop over synthetic replicates
for (m in 1:length(synthetic_list)) {
  syn_data <- synthetic_list[[m]]

  # Loop over attacker models
  for (attacker in attackers) {

    # Loop over tau thresholds
    for (tau in tau_vals) {
      cat("Processing: Replicate", m, "| Model:", attacker, "| τ =", tau, "\n")
      res <- compute_rapid_cat(
        original_data = df,
        synthetic_data = syn_data,
        y = "income",
        tau = tau,
        model = attacker,
        n_boot = n_boot
      )

      results[[length(results) + 1]] <- data.frame(
        Replicate = m,
        Synthesizer = "CART",
        Attacker = toupper(attacker),
        Tau = tau,
        RCIR = res$rcir,
        CI_Lower = res$ci[1],
        CI_Upper = res$ci[2]
      )
    }
  }
}

# Combine results into data.frame
results_df <- do.call(rbind, results)

# Summarize across replicates to build Table 1
library(dplyr)
library(tidyr)

table1 <- results_df %>%
  group_by(Synthesizer, Attacker, Tau) %>%
  summarise(
    RCIR = mean(RCIR),
    CI_Lower = mean(CI_Lower),
    CI_Upper = mean(CI_Upper),
    .groups = "drop"
  ) %>%
  mutate(
    Tau = sprintf("τ = %.2f", Tau),
    Result = sprintf("%.4f (%.4f–%.4f)", RCIR, CI_Lower, CI_Upper)
  ) %>%
  select(Synthesizer, Attacker, Tau, Result) %>%
  pivot_wider(
    names_from = Tau,
    values_from = Result
  ) %>%
  arrange(Synthesizer, factor(Attacker, levels = c("CART", "RF", "GBM", "LOGIT")))

print(table1, row.names = FALSE)
save(table1, result_df, file = "docs/RiskPaper/R/table1_cart.RData")
load("docs/RiskPaper/R/table1_cart.RData")
library(knitr)
library(kableExtra)

# Generate LaTeX table
latex_table <- table1 %>%
  kable(format = "latex", booktabs = TRUE,
        caption = "RAPID-based Relative Confidence Inference Rate (RCIR) by attacker model and threshold.",
        col.names = c("Synthesizer", "Attacker", "$\\tau = 1.10$", "$\\tau = 1.25$", "$\\tau = 1.50$"),
        escape = FALSE, align = "llccc") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

# Print to screen or copy to .tex
cat(latex_table)


library(dplyr)
library(tidyr)

# Step 1: Average RCIR values + CI across replicates
summary_table <- results_df %>%
  group_by(Synthesizer, Attacker, Tau) %>%
  summarise(
    RCIR = mean(RCIR),
    CI_Lower = mean(CI_Lower),
    CI_Upper = mean(CI_Upper),
    .groups = "drop"
  )

# Step 2: Construct main table with formatted results
main_table <- summary_table %>%
  mutate(
    Tau = sprintf("τ = %.2f", Tau),
    Result = sprintf("%.4f (%.4f--%.4f)", RCIR, CI_Lower, CI_Upper)
  ) %>%
  select(Synthesizer, Attacker, Tau, Result) %>%
  pivot_wider(
    names_from = Tau,
    values_from = Result
  )

# Step 3: Compute RCIR max envelope (max RCIR per Synthesizer and Tau)
envelope <- summary_table %>%
  group_by(Synthesizer, Tau) %>%
  summarise(Max_RCIR = max(RCIR), .groups = "drop") %>%
  mutate(
    Tau = sprintf("τ = %.2f", Tau)
  ) %>%
  pivot_wider(
    names_from = Tau,
    values_from = Max_RCIR
  ) %>%
  mutate(
    Attacker = "RCIR max envelope",
    across(starts_with("τ = "), ~ sprintf("\\textbf{%.4f}", .x))
  )

# Step 4: Combine into final table
table1 <- bind_rows(main_table, envelope) %>%
  arrange(Synthesizer, match(Attacker, c("RF", "GBM", "LOGIT", "RCIR max envelope")))


library(knitr)
library(kableExtra)

# Generate LaTeX table
latex_table <- table1 %>%
  kable(format = "latex", booktabs = TRUE,
        caption = "RAPID-based Relative Confidence Inference Rate (RCIR) by attacker model and threshold for the UCI Adult income category (categorical confidential attribute). Values are proportions; 95\\% bootstrap confidence intervals (CIs) in parentheses, computed using 500 non-parametric bootstrap replicates from the original data. For each synthesizer, five synthetic datasets were generated and evaluated. The row \\textit{RCIR max envelope} reports the maximum RCIR point estimate across all attacker models (RF, GBM, LOGIT) for each threshold. This provides a conservative upper bound on attribute-inference risk per synthesizer. \\textcolor{red}{Code is already prepared for extensions to other synthesizer-attacker pairs, e.g., synthpop with ranger, but synthpop with ranger is currently slow.}",
        col.names = c("Synthesizer", "Attacker", "$\\tau = 1.10$", "$\\tau = 1.25$", "$\\tau = 1.50$"),
        escape = FALSE, align = "llccc", longtable = FALSE) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

# Output LaTeX code to console (can redirect to .tex if needed)
cat(latex_table)


#########################
## Now ranger as synth:

# Generate synthetic datasets using ranger
synths_ranger <- synthpop::syn(df, m = 5, method = "ranger")
synthetic_ranger_list <- synths_ranger$syn  # List of 5 synthetic datasets

# Create results container
results_ranger <- list()

# Loop over synthetic replicates
for (m in 1:length(synthetic_ranger_list)) {
  syn_data <- synthetic_ranger_list[[m]]

  # Loop over attacker models
  for (attacker in attackers) {

    # Loop over tau thresholds
    for (tau in tau_vals) {
      cat("Processing: Replicate", m, "| Model:", attacker, "| τ =", tau, "\n")
      res <- compute_rapid_cat(
        original_data = df,
        synthetic_data = syn_data,
        y = "income",
        tau = tau,
        model = attacker,
        n_boot = n_boot
      )

      results_ranger[[length(results_ranger) + 1]] <- data.frame(
        Replicate = m,
        Synthesizer = "CART",
        Attacker = toupper(attacker),
        Tau = tau,
        RCIR = res$rcir,
        CI_Lower = res$ci[1],
        CI_Upper = res$ci[2]
      )
    }
  }
}

# Combine results into data.frame
results_ranger_df <- do.call(rbind, results_ranger)

# Summarize across replicates to build Table 1
library(dplyr)
library(tidyr)

table1_ranger <- results_ranger_df %>%
  group_by(Synthesizer, Attacker, Tau) %>%
  summarise(
    RCIR = mean(RCIR),
    CI_Lower = mean(CI_Lower),
    CI_Upper = mean(CI_Upper),
    .groups = "drop"
  ) %>%
  mutate(
    Tau = sprintf("τ = %.2f", Tau),
    Result = sprintf("%.4f (%.4f–%.4f)", RCIR, CI_Lower, CI_Upper)
  ) %>%
  select(Synthesizer, Attacker, Tau, Result) %>%
  pivot_wider(
    names_from = Tau,
    values_from = Result
  ) %>%
  arrange(Synthesizer, factor(Attacker, levels = c("CART","RF", "GBM", "LOGIT")))

print(table1_ranger, row.names = FALSE)
save(table1_ranger, result_ranger_df,  file = "docs/RiskPaper/R/table1_ranger.RData")

#
# # Optional: Calibration curve of model trained on synthetic data
# attacker_model <- ranger::ranger(
#   formula = income ~ .,
#   data = synthetic_list[[1]],
#   probability = TRUE
# )
#
# # Predict on original data
# probs <- predict(attacker_model, df, type = "response")$predictions
# pred_probs <- probs[, which(colnames(probs) == ">50K")]
# true_labels <- ifelse(df$income == ">50K", 1, 0)
#
# # Calibration plot using caret
# cal <- caret::calibration(as.factor(true_labels) ~ pred_probs, class = "1", cuts = 10)
# autoplot(cal)
# pdf("/Users/matthias/Downloads/calibration_plot.pdf", width = 8, height = 6)
# ggplot2::ggplot(cal)
# dev.off()
