# cmahalanobis.R

utils::globalVariables(c("Species", "Distance", "Comparison"))


#' @importFrom stats mahalanobis
#'
#' @importFrom stats cov
#' @importFrom stats mahalanobis
#' @importFrom stats cov
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom reshape2 melt 
#' @importFrom reshape2 melt
#' 
#' @name cmahalanobis
#' @title Calculate the Mahalanobis distance for each species
#'
#'@description
#'. This function takes a list of data frames as input, where each data frame contains the observations of a species, and returns a matrix with the Mahalanobis distances between each pair of species.
#'
#'
#' @param dataset A list of data frames, where each data frame contains the observations of a species.
#' @param plot Logical, if TRUE, a plot of the Mahalanobis distances is displayed.
#' @param p.value Logical, if TRUE, a matrix of p-values for the distances is returned.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A list containing a matrix with the Mahalanobis distances between each pair of groups, and optionally a matrix of p-values and the plot.
#'
#'
#'
#'
#'
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
#' cmahalanobis(groups, plot = TRUE, p.value = FALSE, plot_title = "Mahalanobis Distance Between Groups")
#' 
#' ceuclide(groups, plot = TRUE, p.value = TRUE, plot_title = "Euclidean Distance Between Groups")
#' cmanhattan(groups, plot = TRUE, p.value = TRUE, plot_title = "Manhattan Distance Between Groups")
#' cchebyshev(groups, plot = TRUE, p.value = TRUE, plot_title = "Chebyshev Distance Between Groups")
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
#' cmahalanobis(groups, plot = TRUE, p.value = TRUE, plot_title = "Mahalanobis Distance Between Groups")
#' ceuclide(groups, plot = TRUE, p.value = TRUE, num.permutations = 1000, plot_title = "Euclidean Distance Between Groups")
#' cmanhattan(groups, plot = TRUE, p.value = TRUE, num.permutations = 1000, plot_title = "Manhattan Distance Between Groups")
#' cchebyshev(groups, plot = TRUE, p.value = TRUE, num.permutations = 1000, plot_title = "Chebyshev Distance Between Groups")
#' 
#'
#' @export
cmahalanobis <- function(dataset, plot = TRUE, p.value = FALSE, plot_title = "Mahalanobis Distance Between Groups") {
  # Verify that the input is a list of data-frames
  if (!is.list(dataset)) {
    stop("The input must be a list of dataframes")
  }
  
  # Verify that each element of the list is a data frame
  for (i in seq_along(dataset)) {
    if (!is.data.frame(dataset[[i]])) {
      stop("Each element of the list must be a dataframe")
    }
  }
  
  # Remove variables not numerical
  dataset <- lapply(dataset, function(df) {
    num_vars <- sapply(df, is.numeric)
    df <- df[, num_vars]
    return(df)
  })
  
  # Replace missing values with arithmetic mean in each dataframe into the list 
  dataset <- lapply(dataset, impute_with_mean)
  
  # Obtain the number of groups
  n <- length(dataset)
  
  # Pre-calculate means and covariances
  means <- lapply(dataset, colMeans)
  covariances <- lapply(dataset, cov)
  
  # Create a empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Mahalanobis distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]]
    cov_i <- covariances[[i]]
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- mean(mahalanobis(x = dataset[[j]], center = mean_i, cov = cov_i))
      }
    }
  }
  
  # Initialize p-value matrix as NULL
  p_values <- NULL
  
  # Calculate p-values only if it is requested by users
  if (p.value) {
    p_values <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      df <- ncol(dataset[[i]]) # Degrees of freedom
      for (j in 1:n) {
        if (i != j) {
          # Calculate p-values using chi-squared distribution
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE)
        }
      }
    }
  }
  
  # If plot is TRUE, call the "plot_mahalanobis_distances" function and print the plot
  if (plot) {
    print(plot_mahalanobis_distances(distances, plot_title))
  }
  
  # Return a list containing distances and optionally p-values
  return(list(distances = distances, p_values = p_values))
}

# Auxiliary function to impute missing values with the mean
impute_with_mean <- function(df) {
  for (i in seq_along(df)) {
    if (any(is.na(df[[i]]) | is.nan(df[[i]]))) {
      df[[i]][is.na(df[[i]]) | is.nan(df[[i]])] <- mean(df[[i]], na.rm = TRUE)
    }
  }
  return(df)
}

# Auxiliary function to print the Mahalanobis distances plot
plot_mahalanobis_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}

#' Calculate Euclidean distance
#'
#' @param dataset A list of dataframes
#' @param plot If TRUE, shows a plot of distances
#' @param p.value If TRUE, calculates p-values for distances
#' @param num.permutations Number of permutations to calculate p-values
#' @param plot_title The title of plot
#' @return A list containing distances and, optionally, p-values
#' @export
ceuclide <- function(dataset, plot = FALSE, p.value = FALSE, num.permutations = 1000, plot_title = "Euclidean Distance Between Groups") {
  # Verify that input is a list of dataframes
  if (!is.list(dataset)) {
    stop("The input must be a list of dataframes")
  }
  
  # Control that each element in the list is a dataframe with at least a numerical variable
  for (df in dataset) {
    if (!is.data.frame(df) || !any(sapply(df, is.numeric))) {
      stop("Ogni elemento della lista deve essere un data frame con almeno una variabile numerica")
    }
  }
  
  # Remove non-numerical variable (Na and NaN) and replace it with arithmetic mean of the variable where they appear 
  dataset <- lapply(dataset, function(df) {
    num_vars <- sapply(df, is.numeric)
    df <- df[, num_vars]
    df <- impute_with_mean(df)
    return(df)
  })
  
  # Pre-calculate mean to each group
  means <- lapply(dataset, colMeans)
  
  # Obtain the number of groups
  n <- length(dataset)
  
  # Create an empty matrix to save distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate euclidean distance between each pair of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Pre-calculated mean of group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- sqrt(sum((mean_i - means[[j]])^2))
      }
    }
  }
  
  # Calculate p-values only if it is requested by users
  if (p.value) {
    p_values <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      mean_i <- means[[i]]  # Pre-calculated mean in the group i
      for (j in 1:n) {
        if (i != j) {
          observed_distance <- distances[i, j]
          permutation_distances <- numeric(num.permutations)
          
          for (k in 1:num.permutations) {
            permuted_data <- sample(unlist(dataset[[j]]))
            permuted_means <- colMeans(matrix(permuted_data, nrow = length(permuted_data)/length(mean_i), byrow = TRUE))
            permutation_distances[k] <- max(abs(mean_i - permuted_means))
          }
          
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        }
      }
    }
  }
  
  # If plot is TRUE, create a plot of distances
  if (plot) {
    print(plot_euclidean_distances(distances, plot_title))
  }
  
  # Return a list containing distances and optionally p-values
  return(list(distances = distances, p_values = p_values))
}

# Auxiliary function to replace Na and NaN values with the mean
impute_with_mean <- function(df) {
  for (i in seq_along(df)) {
    if (is.numeric(df[[i]])) {
      na_index <- is.na(df[[i]]) | is.nan(df[[i]])
      df[[i]][na_index] <- mean(df[[i]], na.rm = TRUE)
    }
  }
  return(df)
}

# Auxiliary function to create a plot of euclidean distances
plot_euclidean_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}


#' Calculate Manhattan distance
#'
#' @param dataset A list of dataframes
#' @param plot If TRUE, show a plot of distances
#' @param p.value If TRUE, calculates p-values to distances
#' @param num.permutations Number of permutations to calculate p-values
#' @param plot_title The title of plot
#' @return A list containing distances and, optionally, p-values
#' @export
# Auxiliary function to create a Manhattan distances plot
cmanhattan <- function(dataset, plot = FALSE, p.value = FALSE, num.permutations = 1000, plot_title = "Manhattan Distance Between Groups") {
  # Verify that input is a list of dataframes
  if (!is.list(dataset)) {
    stop("The input must be a list of dataframes")
  }
  
  # Control that each element in the list is a dataframe with at least a numerical variable
  for (df in dataset) {
    if (!is.data.frame(df) || !any(sapply(df, is.numeric))) {
      stop("Ogni elemento della lista deve essere un data frame con almeno una variabile numerica")
    }
  }
  
  # Remove not numerical variables and replace Na and NaN values with the mean
  dataset <- lapply(dataset, function(df) {
    num_vars <- sapply(df, is.numeric)
    df <- df[, num_vars]
    df <- impute_with_mean(df)
    return(df)
  })
  
  # Pre-calculate the mean to each group
  means <- lapply(dataset, colMeans)
  
  # Obtain the number of groups
  n <- length(dataset)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Manhattan distances between each pair of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Pre-calculated mean in the group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- sum(abs(mean_i - means[[j]]))
      }
    }
  }
  
  # Calculate p-values only if it is requested by users
  if (p.value) {
    p_values <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      mean_i <- means[[i]]  # Pre-calculated mean in the group i
      for (j in 1:n) {
        if (i != j) {
          observed_distance <- distances[i, j]
          permutation_distances <- numeric(num.permutations)
          
          for (k in 1:num.permutations) {
            permuted_data <- sample(unlist(dataset[[j]]))
            permuted_means <- colMeans(matrix(permuted_data, nrow = length(permuted_data)/length(mean_i), byrow = TRUE))
            permutation_distances[k] <- max(abs(mean_i - permuted_means))
          }
          
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        }
      }
    }
  }
  
  # If plot is TRUE, create a plot of distances
  if (plot) {
    print(plot_manhattan_distances(distances, plot_title))
  }
  
  # Return a list containing distances and optionally p-values
  return(list(distances = distances, p_values = p_values))
}

# Auxiliary function to replace Na and NaN values with the mean
impute_with_mean <- function(df) {
  for (i in seq_along(df)) {
    if (is.numeric(df[[i]])) {
      na_index <- is.na(df[[i]]) | is.nan(df[[i]])
      df[[i]][na_index] <- mean(df[[i]], na.rm = TRUE)
    }
  }
  return(df)
}


plot_manhattan_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}



#' Calculate Chebyshev distance
#'
#' @param dataset A list of dataframes
#' @param plot If TRUE, displays a plot of distances
#' @param p.value If TRUE, calculates p-values of distances
#' @param num.permutations Number of permutations to calculate p-values 
#' @param plot_title The title of plot
#' @return A list containing distances and, optionally, p-values
#' @export
cchebyshev <- function(dataset, plot = FALSE, p.value = FALSE, num.permutations = 1000, plot_title = "Chebyshev Distance Between Groups") {
  # Verify that input is a list o dataframes
  if (!is.list(dataset)) {
    stop("The input must be a list of dataframes")
  }
  
  # Control that each element of the list is a dataframes with at least a numerical variable
  for (df in dataset) {
    if (!is.data.frame(df) || !any(sapply(df, is.numeric))) {
      stop("Each element of the list must be dataframes with at least a numerical variable")
    }
  }
  
  # Remove not numerical variables and replace Na and NaN values with the arithmetic mean
  dataset <- lapply(dataset, function(df) {
    num_vars <- sapply(df, is.numeric)
    df <- df[, num_vars]
    df <- impute_with_mean(df)
    return(df)
  })
  
  # Pre-calculate the mean to each group
  means <- lapply(dataset, colMeans)
  
  # Obtain the number of groups
  n <- length(dataset)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Chebyshev distance between each pair of groups
  for (i in 1:n) {
    mean_i <- means[[i]]  # Pre-calculated mean in the group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- max(abs(mean_i - means[[j]]))
      }
    }
  }
  
  # Calculate p-values only if it is requested by users
  if (p.value) {
    p_values <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      mean_i <- means[[i]]  # Pre-calculated mean in the group i
      for (j in 1:n) {
        if (i != j) {
          observed_distance <- distances[i, j]
          permutation_distances <- numeric(num.permutations)
          
          for (k in 1:num.permutations) {
            permuted_data <- sample(unlist(dataset[[j]]))
            permuted_means <- colMeans(matrix(permuted_data, nrow = length(permuted_data)/length(mean_i), byrow = TRUE))
            permutation_distances[k] <- max(abs(mean_i - permuted_means))
          }
          
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        }
      }
    }
  }
  
  # If plot is TRUE, create a plot of distances
  if (plot) {
    print(plot_chebyshev_distances(distances, plot_title))
  }
  
  # Return a list containing distances and optionally p-values
  return(list(distances = distances, p_values = p_values))
}

# Auxiliary function to replace Na and NaN values with the arithmetic mean
impute_with_mean <- function(df) {
  for (i in seq_along(df)) {
    if (is.numeric(df[[i]])) {
      na_index <- is.na(df[[i]]) | is.nan(df[[i]])
      df[[i]][na_index] <- mean(df[[i]], na.rm = TRUE)
    }
  }
  return(df)
}

# Auxiliary function to create a plot of Chebyshev distances
plot_chebyshev_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}


