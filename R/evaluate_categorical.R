#' Evaluate categorical predictions for RAPID risk
#'
#' @description
#' Computes the RAPID metric for a categorical sensitive attribute by comparing
#' predicted class probabilities from a model trained on synthetic data against
#' the true labels in the original data.
#' @param original_data Data frame of original data.
#' @param true_labels A factor vector of true class labels (the sensitive attribute in the original data).
#' @param predicted_probs A matrix or data frame of predicted class probabilities
#'        (rows = observations, columns = classes, typically output from predict()).
#' @param tau A numeric threshold for the relative confidence score (default: 1.0+).
#' @param cat_eval_method Method how the risk score is calculated for categorical variables.
#'
#' @return A list with:
#'   \item{confidence_rate}{Proportion of observations with relative score above \code{tau}.}
#'   \item{relative_score}{Vector of per-observation scores: predicted prob for true class divided by average class baseline.}
#'   \item{true_probs}{Vector of predicted probabilities for the true class label.}
#'
#' @export
evaluate_categorical <- function(true_labels,
                                 predicted_probs,
                                 tau,
                                 original_data,
                                 cat_eval_method = c("RCS", "NCE")) {
  if (!is.factor(true_labels))
    stop("true_labels must be a factor.")
  if (is.data.frame(predicted_probs))
    predicted_probs <- as.matrix(predicted_probs)
  if (is.null(colnames(predicted_probs)))
    stop("predicted_probs must have column names matching the levels of true_labels.")
  if (is.null(cat_eval_method)) {
    cat_eval_method <- "RCS" # default
  } else {
    cat_eval_method <- match.arg(cat_eval_method)
  }

  # g_i: predicted probability for the true class
  g_i <- predicted_probs[cbind(1:nrow(predicted_probs), match(as.character(true_labels), colnames(predicted_probs)))]

  switch(cat_eval_method, RCS = {
    # Option 1 Relative Confidence Score ---------------------------------------
    # observed class levels in both input and predictions
    observed_classes <- intersect(unique(true_labels), colnames(predicted_probs))

    # b_i: baseline probability for each observation based on true class
    b_i <- sapply(unique(as.character(true_labels)), function(k) {
      mean(as.data.frame(predicted_probs)[which(as.character(true_labels) == k), k])
    })[match(as.character(true_labels), observed_classes)]

    # r_i: relative score per observation
    r_i <- g_i / b_i

    # overall confidence rate: proportion of records with relative score > tau
    confidence_rate <- sum(r_i > tau) / length(r_i)

    #  Output reporting table
    result <- data.frame(
      original_data,
      pre_labels = colnames(predicted_probs)[max.col(predicted_probs)] ,
      true_probs = g_i,
      base_prob = b_i,
      relative_score = r_i,
      tau,
      at_risk = ifelse(r_i > tau, TRUE, FALSE)
    )

    result <- result[result$at_risk == TRUE, ]

    # table(colnames(predicted_probs)[max.col(predicted_probs)],
    #       true_labels)
  }, NCE =  {
    # Option 2 Normalized Cross-Entropy ----------------------------------------
    # Cross-entropy (information leaked)
    ce <- -log(pmax(g_i, 1e-10))
    # Normalize by max possible entropy
    n_classes <- ncol(predicted_probs)
    max_entropy <- log(n_classes)
    normalized_ce <- ce / max_entropy
    # Lower CE = better prediction = higher risk
    # Convert to risk score: 1 - normalized_ce
    risk_score <- 1 - normalized_ce
    at_risk <- risk_score > tau

    # overall risk rate: proportion of records being at risk
    risk_rate <- sum(at_risk) / nrow(original_data)

    result <- data.frame(
      original_data,
      predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
      true_prob = g_i,
      risk_score = risk_score,
      at_risk = at_risk
    )
    result <- result[result$at_risk == TRUE, ]
  })

  if (cat_eval_method == "RCS") {
    return(list(
      confidence_rate = confidence_rate,
      relative_score = r_i,
      true_probs = g_i,
      rows_risk_df = result
    ))
  } else {
    return(list(
      risk_rate = risk_rate,
      normalized_ce = normalized_ce,
      true_probs = g_i,
      rows_risk_df = result
    ))
  }
}
