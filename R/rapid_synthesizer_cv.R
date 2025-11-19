rapid_synthesizer_cv <- function(original_data,
                                 synthesizer = NULL,
                                 quasi_identifiers,
                                 sensitive_attribute,
                                 k = 5,
                                 stratified = TRUE,
                                 return_details = FALSE,
                                 return_all_records = FALSE,
                                 seed = 2025,
                                 trace = TRUE,
                                 ...) {
  #  Setup
  if (!is.null(seed))
    set.seed(seed)

  if (is.null(synthesizer)) {
    if (!requireNamespace("synthpop", quietly = TRUE)) {
      stop("synthpop required or provide custom synthesizer")
    }
    synthesizer <- function(data, seed = NULL) {
      synthpop::syn(data,
                    m = 1,
                    seed = seed,
                    print.flag = FALSE)$syn
    }
  }

  stopifnot(k >= 2, nrow(original_data) >= k, is.function(synthesizer))
  target_vec <- original_data[[sensitive_attribute]]

  is_categorical <- is.factor(target_vec) ||
    is.character(target_vec)

  # Create Folds
  if (stratified &&
      is_categorical && requireNamespace("caret", quietly = TRUE)) {
    folds <- create_stratified_folds(target_vec, k) # utils function
  } else {
    folds <- split(1:nrow(original_data), sample(rep(1:k, length.out = nrow(original_data))))
  }

  #  Run CV
  if (trace) {
    cat(sprintf("  RAPID Synthesizer CV (%d folds", k))
    if (stratified && is_categorical)
      cat(", stratified")
  }

  # Pre-allocate results
  fold_results <- vector("list", k)
  fold_data <- if (return_all_records)
    vector("list", k)
  else
    NULL

  for (fold_id in seq_len(k)) {
    if (trace)
      cat(sprintf("  Fold %d/%d: ", fold_id, k))

    # Split
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(original_data)), test_idx)

    # Synthesize on train
    if (trace)
      cat(sprintf("synthesizing %d → ", length(train_idx)))

    if (!is.null(seed)) {
      synth <- synthesizer(original_data[train_idx, ], seed = seed + fold_id)
    } else {
      synth <- synthesizer(original_data[train_idx, ])
    }

    # Evaluate on test
    if (trace)
      cat(sprintf("evaluating %d... ", length(test_idx)))


    result <- rapid(
      original_data = original_data[test_idx, ],
      synthetic_data = synth,
      quasi_identifiers = quasi_identifiers,
      sensitive_attribute = sensitive_attribute,
      seed = NULL,
      trace = FALSE,
      return_all_records = return_all_records,
      ...
    )

    # Extract confidence rate
    conf_rate <- result$risk$confidence_rate

    if (trace)
      cat(sprintf("✓ (risk: %.3f)\n", conf_rate))

    # Store results as data.frame directly
    if (is_categorical) {
      fold_results[[fold_id]] <- data.frame(
        fold = fold_id,
        confidence_rate = conf_rate,
        accuracy = result$metrics$accuracy,
        stringsAsFactors = FALSE
      )
    } else {
      fold_results[[fold_id]] <- data.frame(
        fold = fold_id,
        confidence_rate = conf_rate,
        mae = result$metrics$mae,
        rmse = result$metrics$rmse,
        rmae = result$metrics$rmae,
        rrmse = result$metrics$rrmse,
        stringsAsFactors = FALSE
      )
    }


    if (return_all_records) {
      fold_data[[fold_id]] <- result$risk$rows_risk_df
      fold_data[[fold_id]]$fold <- fold_id
    }
  }
  # ===== Aggregate =====
  fold_df <- do.call(rbind, fold_results)
  conf_rates <- fold_df$confidence_rate

  cv_summary <- list(
    mean = mean(conf_rates),
    sd = sd(conf_rates),
    se = sd(conf_rates) / sqrt(k),
    median = median(conf_rates),
    min = min(conf_rates),
    max = max(conf_rates),
    ci_lower = mean(conf_rates) - qnorm(0.975) * sd(conf_rates) / sqrt(k),
    ci_upper = mean(conf_rates) + qnorm(0.975) * sd(conf_rates) / sqrt(k)
  )

  if (is_categorical) {
    cv_summary$mean_accuracy <- mean(fold_df$accuracy, na.rm = TRUE)
    cv_summary$sd_accuracy <- sd(fold_df$accuracy, na.rm = TRUE)
  } else {
    cv_summary$mean_mae <- mean(fold_df$mae, na.rm = TRUE)
    cv_summary$mean_rmse <- mean(fold_df$rmse, na.rm = TRUE)
    cv_summary$mean_rmae <- mean(fold_df$rmae, na.rm = TRUE)
    cv_summary$mean_rrmse <- mean(fold_df$rrmse, na.rm = TRUE)
  }

  result <- list(
    cv_summary = cv_summary,
    cv_details = if (return_details)
      fold_df
    else
      NULL,
    fold_data = if (return_all_records)
      do.call(rbind, fold_data)
    else
      NULL,
    settings = list(
      k = k,
      cv_type = "synthesizer",
      stratified = stratified,
      is_categorical = is_categorical,
      n_original = nrow(original_data),
      sensitive_attribute = sensitive_attribute
    )
  )

  class(result) <- c("rapid_synthesizer_cv", "rapid_cv", "list")

  if (trace) {
    cat(
      sprintf(
        "  Mean Risk: %.4f [%.4f, %.4f]\n",
        cv_summary$mean,
        cv_summary$ci_lower,
        cv_summary$ci_upper
      )
    )
  }

  return(result)
}
