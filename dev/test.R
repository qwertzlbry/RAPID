############################################################
## TESTSKRIPT für rapid()
## Originaldaten + synthetische Daten mit synthpop
############################################################

## 1) Pakete
# install.packages(c("synthpop", "ranger", "gbm", "rpart"))
library(synthpop)
devtools::load_all()

## 2) Beispieldaten laden
data("SD2011", package = "synthpop")

## Variablenauswahl (geeignet für rapid)
vars <- c("sex", "age", "edu", "marital", "income", "ls", "wkabint")
orig_df <- SD2011[, vars]

## Missing-Codierung bereinigen
orig_df$income[orig_df$income == -8] <- NA

## 3) Synthetische Daten erzeugen
set.seed(123)

syn_obj <- syn(orig_df)
syn_df  <- syn_obj$syn

## 4) Definition von QID & Sensitive Variable
quasi_identifiers   <- c("sex", "age", "edu", "marital")
sensitive_attribute <- "income"   # numerisch → nutzt numeric-Teil

## 5) rapid() laufen lassen (NUMERISCHER TEST)
risk_num <- rapid(
  original_data      = orig_df,
  synthetic_data     = syn_df,
  quasi_identifiers  = quasi_identifiers,
  sensitive_attribute = sensitive_attribute,

  model_type = "rf",

  num_na_strategy    = "median",
  num_error_metric   = "rmse",
  num_epsilon_type   = "Percentage",
  num_epsilon        = 0.1,

  seed  = 2025,
  trace = TRUE
)

risk_num$risk

sensitive_attribute_cat <- "marital"
risk_cat <- rapid(
  original_data      = orig_df,
  synthetic_data     = syn_df,
  quasi_identifiers  = c("sex", "age", "edu"),
  sensitive_attribute = sensitive_attribute_cat,

  model_type = "gbm",   # passend für kategorisch

  cat_eval_method = "RCS_marginal",
  cat_tau         = 0.5,

  seed  = 2025,
  trace = TRUE
)


risk_cat$risk$confidence_rate


models_num <- c("lm", "rf", "gbm", "cart")

rows_all_num <- do.call(rbind, lapply(models_num, function(m) {

  res <- rapid(
    original_data       = orig_df,
    synthetic_data      = syn_df,
    quasi_identifiers   = c("sex", "age", "edu"),
    sensitive_attribute = "income",
    model_type          = m,
    num_epsilon         = 0.1,
    seed                = 2025
  )

  df <- res$risk$rows_risk_df
  df$model <- m
  df
}))









models_num <- c("lm", "rf", "gbm", "cart")

rows_all_num <- do.call(rbind, lapply(models_num, function(m) {

  res <- rapid(
    original_data       = orig_df,
    synthetic_data      = syn_df,
    quasi_identifiers   = c("sex", "age", "edu"),
    sensitive_attribute = "income",
    model_type          = m,

    # ✅ WICHTIG: für numerisch zwingend notwendig
    num_epsilon_type = "Percentage",   # oder "Value"
    num_epsilon      = 0.1,             # 10% Toleranz (oder z.B. 5 bei "Value")

    num_error_metric = "rmse",          # optional aber sinnvoll
    seed             = 2025
  )

  df <- res$risk$rows_risk_df
  df$model <- m
  df
}))

ggplot(rows_all_num, aes(x = model, y = relative_error_percent)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_minimal() +
  labs(
    title = "Numeric Disclosure Risk – Absolute Error",
    x = "Model",
    y = "Absolute Error"
  )

## KAT
models_cat <- c("rf", "gbm", "cart")   # logit bewusst NICHT hier

comparison_cat <- do.call(rbind, lapply(models_cat, function(m) {

  res <- rapid(
    original_data       = orig_df,
    synthetic_data      = syn_df,
    quasi_identifiers   = c("sex", "age", "edu"),
    sensitive_attribute = "marital",

    model_type      = m,
    cat_eval_method = "RCS_marginal",
    cat_tau         = 0.5,
    seed            = 2025
  )

  extract_rapid_metrics(
    res$risk,
    model       = m,
    sensitive   = "marital",
    eval_method = "RCS_marginal",
    type        = "categorical"
  )
}))


library(ggplot2)

ggplot(comparison_cat, aes(x = model, y = normalized_gain)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "RAPID Model Comparison (Categorical)",
    y = "Confidence Rate",
    x = "Model"
  )


