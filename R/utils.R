create_stratified_folds <- function(y, k) {
  folds <- vector("list", k)

  # For each class/level
  for (class in unique(y)) {
    class_indices <- which(y == class)
    n_class <- length(class_indices)

    # Shuffle and assign to folds
    shuffled <- sample(class_indices)
    fold_assignment <- rep(1:k, length.out = n_class)

    # Add to folds
    for (fold_id in 1:k) {
      folds[[fold_id]] <- c(folds[[fold_id]], shuffled[fold_assignment == fold_id])
    }
  }

  # Shuffle within each fold (optional, for randomness)
  folds <- lapply(folds, sample)

  return(folds)
}


# y <- factor(c(rep("A", 100), rep("B", 200), rep("C", 50)))
# folds <- create_stratified_folds(y, k = 5)
# lapply(folds, function(idx) table(y[idx]))

