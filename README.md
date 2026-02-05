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
devtools::install_github("yourusername/rapid")
```

## Quick Start
```r
library(rapid)

# Assess risk for a categorical sensitive attribute
result <- rapid(
  original_data = original_df,
  synthetic_data = synthetic_df, 
  sensitive_attribute = "income",
  quasi_identifiers = c("age", "education", "occupation"),
  cat_tau = 0.3  
)

# View summary and plot
summary(result)
plot(result)
```

## Citation

If you use RAPID in your research, please cite:
```
Templ, M., Thees, O., & Müller, R. (2025). 
RAPID: Risk of Attribute Prediction-Induced Disclosure in Synthetic Microdata.
arXiv preprint [TBA]
```

## Documentation

- **Paper**: [arXiv link - TBA]
- **Documentation**: [TBA]
- **Issues**: [TBA]
