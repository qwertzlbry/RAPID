# RAPID: Risk Assessment through Prediction for Inference-based Disclosure

RAPID is an R package for assessing attribute inference risk in synthetic microdata using machine learning-based attacker models. 
It provides record-level risk scores with interpretable thresholds for both categorical and continuous sensitive attributes.

## Key Features
- **Model-based risk assessment**: Simulates sophisticated attackers using machine learning methods
- **Unified framework**: Handles both categorical (normalized gain) and continuous (symmetric relative error) sensitive attributes
- **Baseline-corrected scoring**: Accounts for class imbalance and marginal predictability
- **Record-level risk flags**: Identifies which individuals are at highest inference risk
- **Mixed-type data support**: Works naturally with continuous, categorical, and mixed quasi-identifiers
- **Diagnostic tools**: Attribution analysis to identify which quasi-identifier patterns drive risk
- **Flexible evaluation**: Supports both match-based and holdout-based assessment modes


## Installation
```r
# Install from GitHub
devtools::install_github("qwertzlbry/RAPID")
```

## Quick Start
```r
library(RAPID)
library(synthpop)

# Simulate original data with dependencies
set.seed(2025)
n <- 1000

age <- sample(20:70, n, replace = TRUE)
education <- factor(sample(c("low", "medium", "high"), n, replace = TRUE))
gender <- factor(sample(c("M", "F"), n, replace = TRUE))

# Disease status depends on age and education
disease_status <- sapply(1:n, function(i) {
  probs <- c(0.6, 0.2, 0.2)  # Base: healthy, diabetic, hypertensive
  if (age[i] > 55) probs <- c(0.3, 0.5, 0.2)
  if (age[i] > 55 && education[i] == "low") probs <- c(0.2, 0.6, 0.2)
  sample(c("healthy", "diabetic", "hypertensive"), 1, prob = probs)
})

data_orig <- data.frame(age, education, gender, 
                        disease_status = factor(disease_status))

# Generate synthetic data
data_syn <- syn(data_orig, method = "cart", seed = 2025)$syn


# Assess risk for a categorical sensitive attribute
result <- rapid(
  original_data = data_orig,
  synthetic_data = data_syn,
  quasi_identifiers = c("age", "education", "gender"),
  sensitive_attribute = "disease_status",
  model_type = "rf",
  cat_tau = 0.3
)

# View print, summary and plot
print(result)
summary(result)
plot(result)
```

## Citation

If you use RAPID in your research, please cite:
```
Templ, M., Thees, O., & Müller, R. (2026). 
RAPID: Risk of Attribute Prediction-Induced Disclosure in Synthetic Microdata.
arXiv preprint [TBA]
```

## Documentation

- **Paper**: [arXiv link - TBA]
- **Documentation**: [TBA]
- **Issues**: [TBA]
