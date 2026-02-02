#' @title Evaluate Categorical Attribute Inference Risk (RAPID)
#'
#' @description
#' Computes the RAPID metric for a categorical sensitive attribute by comparing
#' predicted class probabilities from a model trained on synthetic data against
#' the true labels in the original data. Supports multiple evaluation methods
#' for measuring inferential disclosure risk.
#'
#' @param A A factor vector of true class labels (the sensitive attribute in the original data).
#' @param B A matrix or data frame of predicted class probabilities
#'        (rows = observations, columns = classes, typically output from predict()).
#' @param cat_tau A numeric threshold for the risk score. Interpretation depends on method:
#'   \itemize{
#'     \item For \code{RCS_conditional}: threshold for ratio (typically 1.25)
#'     \item For \code{RCS_marginal}: threshold for normalized gain (typically 0.3)
#'     \item For \code{NCE}: threshold for risk score (typically 0.5-0.7)
#'   }
#' @param original_data Data frame of original data.
#' @param cat_eval_method Method for calculating risk score. Options:
#'   \describe{
#'     \item{\code{RCS_conditional}}{Relative Confidence Score with class-conditional baseline.
#'           Measures if observation is an outlier within its class. Uses simple ratio \eqn{r_i = g_i / b_i}
#'           where \eqn{b_i} is the mean prediction for observations in the same true class.}
#'     \item{\code{RCS_marginal}}{Relative Confidence Score with marginal baseline (recommended).
#'           Measures if attribute can be inferred better than baseline rate from original data.
#'           Uses normalized gain \eqn{(g_i - b_i) / (1 - b_i)} where \eqn{b_i} is the marginal
#'           frequency of the true class in original data. Provides fair comparison across classes.}
#'     \item{\code{NCE}}{Normalized Cross-Entropy. Measures information leakage using
#'           entropy-based approach.}
#'   }
#'   Default: \code{"RCS_marginal"}
#' @param sensitive_attribute Character string specifying the name of the sensitive
#'        attribute column in original_data (required for RCS_marginal).
#' @param return_all_records Logical; if \code{TRUE}, return all records with their risk
#'   status. If \code{FALSE} (default), return only at-risk records.
#'
#' @return A list with elements depending on the method. All methods return:
#' \describe{
#'   \item{method}{Character string indicating which method was used}
#'   \item{confidence_rate (or risk_rate for NCE)}{Proportion of observations at risk (0-1)}
#'   \item{n_at_risk}{Number of records flagged as at-risk}
#'   \item{percentage}{Percentage of at-risk records (0-100)}
#'   \item{rows_risk_df}{Data frame with all records (if \code{return_all_records = TRUE})
#'     or only at-risk records (if \code{FALSE}), including original data, predictions,
#'     risk scores, and \code{at_risk} flag}
#' }
#'
#'   Method-specific additional components:
#'
#'   For \code{RCS_conditional}:
#'   \itemize{
#'     \item{\code{relative_score}}{Vector of per-observation ratios: \eqn{g_i / b_i}}
#'     \item{\code{baseline}}{Vector of class-conditional baselines}
#'     \item{\code{true_probs}}{Vector of predicted probabilities for the true class (\eqn{g_i})}
#'   }
#'
#'   For \code{RCS_marginal}:
#'   \itemize{
#'     \item{\code{normalized_gain}}{Vector of per-observation normalized improvements over baseline}
#'     \item{\code{baseline}}{Vector of marginal baselines from original data}
#'     \item{\code{true_probs}}{Vector of predicted probabilities for the true class (\eqn{g_i})}
#'   }
#'
#'   For \code{NCE}:
#'   \itemize{
#'     \item{\code{normalized_ce}}{Vector of normalized cross-entropy values}
#'     \item{\code{true_probs}}{Vector of predicted probabilities for the true class}
#'   }
#'
#' @details
#' The choice of method affects both the interpretation and recommended threshold:
#'
#' \strong{RCS_conditional} (class-conditional baseline):
#' \itemize{
#'   \item Measures: "Is this observation unusual compared to others in its class?"
#'   \item Use case: Identifying outliers within classes
#'   \item Recommended cat_tau: 1.25 (25\% better than class average)
#'   \item Note: Risk may decrease as model quality increases (paradoxical behavior)
#' }
#'
#' \strong{RCS_marginal} (marginal baseline, recommended):
#' \itemize{
#'   \item Measures: "Can the attribute be inferred better than baseline rate?"
#'   \item Use case: Inferential disclosure risk assessment
#'   \item Recommended cat_tau: 0.3 (30\% of possible improvement over baseline)
#'   \item Advantages: Fair across classes, accounts for class imbalance, consistent with
#'     numeric case, monotonic with model quality
#' }
#'
#' \strong{NCE} (Normalized Cross-Entropy):
#' \itemize{
#'   \item Measures: Information leakage via entropy
#'   \item Use case: Information-theoretic risk assessment
#'   \item Recommended cat_tau: 0.5-0.7
#' }
#'
#' @examples
#' \dontrun{
#' # Example with RCS_marginal (recommended)
#' library(ranger)
#'
#' # Fit model on synthetic data
#' fit <- ranger(disease_status ~ age + income + education,
#'               data = synthetic_data, probability = TRUE)
#'
#' # Get predictions on original data
#' pred_probs <- predict(fit, original_data)$predictions
#'
#' # Evaluate risk
#' result <- evaluate_categorical(
#'   A = original_data$disease_status,
#'   B = pred_probs,
#'   cat_tau = 0.3,
#'   original_data = original_data,
#'   sensitive_attribute = "disease_status",
#'   cat_eval_method = "RCS_marginal",
#'   return_all_records = TRUE
#' )
#'
#' # Access results
#' print(result$confidence_rate)  # Proportion at risk
#' print(result$n_at_risk)        # Number of at-risk records
#' head(result$rows_risk_df)      # Per-record details
#'
#' # Plot normalized gain distribution
#' hist(result$normalized_gain,
#'      main = "Distribution of Normalized Gain",
#'      xlab = "Normalized Gain")
#' abline(v = 0.3, col = "red", lty = 2)  # Threshold
#'
#' # Example with RCS_conditional
#' result_cond <- evaluate_categorical(
#'   A = original_data$disease_status,
#'   B = pred_probs,
#'   cat_tau = 1.25,
#'   original_data = original_data,
#'   sensitive_attribute = "disease_status",
#'   cat_eval_method = "RCS_conditional"
#' )
#'
#' # Example with NCE
#' result_nce <- evaluate_categorical(
#'   A = original_data$disease_status,
#'   B = pred_probs,
#'   cat_tau = 0.6,
#'   original_data = original_data,
#'   sensitive_attribute = "disease_status",
#'   cat_eval_method = "NCE"
#' )
#' }
#'
#' @export
evaluate_categorical <- function(A,
                                 B,
                                 cat_tau,
                                 original_data,
                                 cat_eval_method = c("RCS_conditional", "RCS_marginal", "NCE"),
                                 sensitive_attribute,
                                 return_all_records = FALSE) {
  true_labels <- A
  predicted_probs <- B

  if (!is.factor(true_labels))
    stop("true_labels must be a factor.")
  if (is.data.frame(predicted_probs))
    predicted_probs <- as.matrix(predicted_probs)
  if (is.null(colnames(predicted_probs)))
    stop("predicted_probs must have column names matching the levels of true_labels.")
  if (is.null(cat_eval_method)) {
    cat_eval_method <- "RCS_marginal" # default
  } else {
    cat_eval_method <- match.arg(cat_eval_method)
  }

  # g_i: predicted probability for the true class
  g_i <- predicted_probs[cbind(1:nrow(predicted_probs), match(as.character(true_labels), colnames(predicted_probs)))]

    switch(cat_eval_method,

    RCS_conditional = {
                    # Option 1 Relative Confidence Score ---------------------------------------
                    # observed class levels in both input and predictions
                    observed_classes <- intersect(unique(true_labels), colnames(predicted_probs))

                    # Original ---
                    # b_i: baseline probability for each observation based on true class
                    # b_i <- sapply(unique(as.character(true_labels)), function(k) {
                    #   mean(as.data.frame(predicted_probs)[which(as.character(true_labels) == k), k])
                    # })[match(as.character(true_labels), observed_classes)]

                    # FIX 1 ---
                    # Compute baseline for each class
                    baseline_by_class <- sapply(unique(as.character(true_labels)), function(k) {
                      mean(predicted_probs[which(as.character(true_labels) == k), k])
                    })
                    names(baseline_by_class) <- unique(as.character(true_labels))

                    # Map baseline to each observation
                    b_i <- baseline_by_class[as.character(true_labels)]

                    #FIX 2 Marginal Baseline to be implemented ---
                    # n_classes <- ncol(predicted_probs)
                    # b_i <- rep(1 / n_classes, nrow(predicted_probs))


                    # r_i: relative score per observation
                    r_i <- g_i / b_i

                    # overall confidence rate: proportion of records with relative score > cat_tau
                    confidence_rate <- sum(r_i > cat_tau) / length(r_i)

                    #  Output reporting table
                    result <- data.frame(
                      original_data,
                      pre_labels = colnames(predicted_probs)[max.col(predicted_probs)] ,
                      true_probs = g_i,
                      base_prob = b_i,
                      relative_score = r_i,
                      cat_tau,
                      at_risk = ifelse(r_i > cat_tau, TRUE, FALSE),
                      stringsAsFactors = FALSE
                    )
                    if (return_all_records) {
                      result
                    } else {
                      result <- result[as.logical(result$at_risk == TRUE), ]
                    }

                  # table(colnames(predicted_probs)[max.col(predicted_probs)],
                  #       true_labels)
  }, RCS_marginal = {
                  # Marginal baseline from ORIGINAL dataset
                  marginal_freq_original <- prop.table(table(original_data[[sensitive_attribute]]))
                  b_i <- as.numeric(marginal_freq_original[as.character(true_labels)])

                  # Normalized improvement (0 to 1)
                  max_improvement <- 1 - b_i
                  normalized_gain <- (g_i - b_i) / max_improvement

                  # At risk when gain > cat_tau
                  at_risk <- normalized_gain > cat_tau

                  confidence_rate <- sum(at_risk) / length(at_risk)

                  result <- data.frame(
                    original_data,
                    predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
                    true_prob = g_i,
                    baseline = b_i,
                    normalized_gain = normalized_gain,
                    cat_tau = cat_tau,
                    at_risk = at_risk,
                    stringsAsFactors = FALSE
                  )

                  if (return_all_records) {
                    result
                  } else {
                    result <- result[as.logical(result$at_risk == TRUE), ]
                  }


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
                  at_risk <- risk_score > cat_tau

                  # overall risk rate: proportion of records being at risk
                  risk_rate <- sum(at_risk) / nrow(original_data)

                  result <- data.frame(
                    original_data,
                    predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
                    true_prob = g_i,
                    risk_score = risk_score,
                    at_risk = at_risk,
                    stringsAsFactors = FALSE
                  )
                  if (return_all_records) {
                    result
                  } else {
                    result <- result[as.logical(result$at_risk == TRUE), ]
                  }
                })

  if (cat_eval_method == "RCS_conditional") {
    return(list(
      method = "RCS_conditional",
      confidence_rate = confidence_rate,
      n_at_risk = sum(result$at_risk),
      percentage = confidence_rate * 100,
      relative_score = r_i,
      baseline = b_i,
      true_probs = g_i,
      rows_risk_df = result
    ))
  } else if (cat_eval_method == "RCS_marginal") {
    return(list(
      method = "RCS_marginal",
      confidence_rate = confidence_rate,
      n_at_risk = sum(result$at_risk),
      percentage = confidence_rate * 100,
      normalized_gain = as.numeric(normalized_gain),
      baseline = b_i,
      true_probs = g_i,
      rows_risk_df = result
    ))
  } else {  # NCE
    return(list(
      method = "NCE",
      risk_rate = risk_rate,
      normalized_ce = normalized_ce,
      true_probs = g_i,
      rows_risk_df = result
    ))
  }
}
