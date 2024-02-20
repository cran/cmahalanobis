# cmahalanobis.R

#' @importFrom stats mahalanobis
#'
#' @importFrom stats cov
#'
#' @name cmahalanobis
#' @title Calculate the Mahalanobis distance for each species
#'
#'@description
#'. This function takes a list of data frames as input, where each data frame contains the observations of a species, and returns a matrix with the Mahalanobis distances between each pair of species.
#'
#'
#' @param dataset A list of data frames, where each data frame contains the observations of a species.
#'
#' @return A matrix with the Mahalanobis distances between each pair of species.
#' @examples
#' # Example with the iris dataset
#' library(stats)
#' # Split the data into 3 parts for each species
#' setosa <- subset(iris, Species == "setosa")
#' setosa <- setosa[,-5]
#' versicolor <- subset(iris, Species == "versicolor")
#' versicolor <- versicolor[,-5]
#' virginica <- subset(iris, Species == "virginica")
#' virginica <- virginica[,-5]
#'
#' # Create a list with the three groups of flowers
#' groups <- list(setosa, versicolor, virginica)
#'
#' # Calculate the Mahalanobis distance with the cmahalanobis function
#' cmahalanobis(groups)
#'
#' # Example with the mtcars dataset
#' library(stats)
#' # Split the data into 2 parts for each type of transmission
#' auto <- subset(mtcars, am == 0)
#' auto <- auto[,-9]
#' manual <- subset(mtcars, am == 1)
#' manual <- manual[,-9]
#'
#' # Create a list with the two groups of cars
#' groups <- list(auto, manual)
#'
#' # Calculate the Mahalanobis distance with the cmahalanobis function
#' cmahalanobis(groups)
#'
#' @export
cmahalanobis <- function(dataset) {
  # Check that the input is a list
  if (!is.list(dataset)) {
    stop("The input must be a list of data frames")
  }

  # Check that each element of the list is a data frame
  for (i in seq_along(dataset)) {
    if (!is.data.frame(dataset[[i]])) {
      stop("Each element of the list must be a data frame")
    }
  }

  # Remove the non-numeric variables from each data frame in the list
  dataset <- lapply(dataset, function(df) {
    # Find the numeric variables in the data frame
    num_vars <- sapply(df, is.numeric)

    # Keep only the numeric variables in the data frame
    df <- df[, num_vars]

    # Return the updated data frame
    return(df)
  })

  # Replace the missing values with 0 in each data frame in the list
  dataset <- lapply(dataset, function(df) {
    # Find the missing values in the data frame
    na_vals <- is.na(df)

    # Replace the missing values with 0 in the data frame
    df[na_vals] <- 0

    # Return the updated data frame
    return(df)
  })

  # Get the number of groups
  n <- length(dataset)

  # Create an empty matrix to save the distances
  distances <- matrix(0, nrow = n, ncol = n)

  # Calculate the Mahalanobis distance between each pair of groups
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the mean vector and the covariance matrix of group i
      mean_i <- colMeans(dataset[[i]])
      cov_i <- cov(dataset[[i]])

      # Calculate the Mahalanobis distance between group i and group j
      distances[i, j] <- mean(mahalanobis(x = dataset[[j]], center = mean_i, cov = cov_i))
    }
  }

  # Return the matrix of distances
  return(distances)
}

