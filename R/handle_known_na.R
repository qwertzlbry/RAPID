#' @title Handle Missing Values in Quasi-Identifiers
#'
#' @description
#' For each quasi-identifier, this function handles missing values by:
#' - For factors: converting \code{NA} to a new level \code{"missing"}.
#' - For numeric variables: stopping with an error if \code{NA} values are present.
#'
#' This function ensures that modeling can proceed without errors due to missing values in quasi-identifiers.
#'
#' @param data A \code{data.frame} containing the data to check.
#' @param quasi_identifiers A character vector specifying the quasi-identifier variables.
#'
#' @return A modified \code{data.frame} with missing values handled.
#' @export
handle_known_na <- function(data, quasi_identifiers) {
  for (var in quasi_identifiers) {
    if (!var %in% names(data)) next
    
    if (is.factor(data[[var]])) {
      if (anyNA(data[[var]])) {
        data[[var]] <- addNA(data[[var]])
        levels(data[[var]])[is.na(levels(data[[var]]))] <- "missing"
        warning(sprintf("NA in factor '%s' converted to level 'missing'.", var))
      }
    } else if (is.numeric(data[[var]]) && anyNA(data[[var]])) {
      stop(sprintf("Variable '%s' is numeric and contains NA. Please handle missing values before modeling.", var))
    }
  }
  
  return(data)
}