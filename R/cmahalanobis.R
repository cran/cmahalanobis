# cmahalanobis.R

utils::globalVariables(c("Species", "Distance", "Comparison"))
utils::globalVariables(c("Var1", "Var2", "value", "element_text"))
utils::globalVariables(c("p_values", "distances"))

#' @importFrom stats mahalanobis
#' @importFrom mice mice
#' @importFrom mice complete
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
#' @title Calculate the Mahalanobis distance for each species.
#'
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Mahalanobis distances about it.
#'
#'
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Mahalanobis distances matrix.
#' @param plot Logical, if TRUE, a plot of Mahalanobis distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Mahalanobis distances between each pair of groups and the plot.
#'
#'
#' @examples
#' # Example with the iris dataset
#'
#' data(iris)
#' 
#' # Calculate the Mahalanobis distance with the cmahalanobis function
#' cmahalanobis(iris, ~Species, plot = TRUE, plot_title = "Mahalanobis Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' data(mtcars)
#' 
#' # Calculate the Mahalanobis distance with the cmahalanobis function
#' cmahalanobis(mtcars, ~am, plot = TRUE, plot_title = "Mahalanobis Distance Between Groups")
#' 
#' @export
cmahalanobis <- function(dataset, formula, plot = TRUE, plot_title = "Mahalanobis Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variables in each group
  groups <- lapply(groups, function(df) {
    df <- df[, sapply(df, is.numeric)]  # Select only numeric columns
    return(df)
  })
  
  # Replace missing values with arithmetic mean in each dataframe into the list
  groups <- lapply(groups, impute_with_mean)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means and covariances
  means <- lapply(groups, colMeans)
  covariances <- lapply(groups, function(df) {
    cov_matrix <- cov(df)
    n <- nrow(cov_matrix)
    reg <- 0.01  # Regularization value
    cov_matrix <- cov_matrix + diag(reg, nrow = n)
    return(cov_matrix)
  })
  
  # Create a empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Mahalanobis distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]]
    cov_i <- covariances[[i]]
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- mean(mahalanobis(x = groups[[j]], center = mean_i, cov = cov_i))
      }
    }
  }
  
  
  # If plot is TRUE, call the "plot_mahalanobis_distances" function and print the plot
  if (plot) {
    print(plot_mahalanobis_distances(distances, plot_title))
  }
  
  # Return a list containing distances 
  return(list(distances = distances))
}

impute_with_mean <- function(df) {
  # Calculate mean for each columns, ignoring missing values
  means <- colMeans(df, na.rm = TRUE)
  
  # Replace missing values with the corresponding mean
  for (i in 1:ncol(df)) {
    df[is.na(df[, i]), i] <- means[i]
  }
  
  return(df)
}

# Auxiliary function to print the Mahalanobis distances plot 
plot_mahalanobis_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- reshape2::melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}






#' @name generate_report_cmahalanobis
#' @title Generate a Microsoft Word document about Mahalanobis distance matrix and p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about Mahalanobis distance matrix and p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Mahalanobis distances matrix and p_values matrix.
#' @param pvalue.method A method with which you want to calculate pvalue matrix.The default method is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations A number of permutations to define if you choose "permutation".
#' @param num.bootstraps A number of bootstrap to define if you choose "bootstrap".
#' @return A Microsoft Word document about Mahalanobis distances matrix and p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cmahalanobis(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cmahalanobis(mtcars, ~am)
#' 
#' @export
generate_report_cmahalanobis <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cmahalanobis_results <- cmahalanobis(dataset, formula)
  distances <- cmahalanobis_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "chisq")  # Adjust method if needed
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "permutation")  # Adjust method if needed
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "bootstrap")  # Adjust method if needed
  }
  
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "reportcmahalanobis.docx")
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Assuming the path to your template is correct
  rmarkdown::render(system.file("rmarkdown", "template_report_cmahalanobis.Rmd", package = "cmahalanobis"),
                    params = list(distances = distances, p_values = p_values),
                    output_file = "reportcmahalanobis.docx")
}



#' @name pvaluescmaha
#' @title Calculate p_values matrix for each species, using Mahalanobis distance as a base.
#'
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of the species and a plot if the user select TRUE using Mahalanobis distance for distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Mahalanobis distances matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq".Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot a p_values heatmap. The default value is TRUE.
#' @return A list containing the p-values matrix and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescmaha(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescmaha(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescmaha <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
  # Verify that the input is a dataframe
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the name of response variable by the formula
  response <- all.vars(formula)[1]
  
  # Verify that the response variable exists in the dataset
  if (!response %in% names(dataset)) {
    stop(paste("Response variable", response, "is absent in the dataset"))
  }
  
  # Calculate Mahalanobis distances using "cmahalanobis" function
  mahalanobis_results <- cmahalanobis(dataset, formula, plot = FALSE)
  distances <- mahalanobis_results$distances
  
  # Obtain groups number
  n <- nrow(distances)
  
  # Initialize p_values matrix
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate p_values
  for (i in 1:n) {
    df <- ncol(dataset) - 1  # Degrees of freedom (adjusted for covariance matrix estimate)
    for (j in 1:n) {
      if (i != j) {
        # Choose p_values method calculation based user input
        if (pvalue.method == "chisq") {
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE, log.p = TRUE)
        } else if (pvalue.method == "permutation") {
          # Permutation test
          observed_distance <- distances[i, j]
          permutation_distances <- replicate(num.permutations, {
            # Permute labels group and recalculate the distance
            permuted_data <- dataset
            permuted_data[[response]] <- sample(dataset[[response]])
            permuted_results <- cmahalanobis(permuted_data, formula, plot = FALSE)
            permuted_distances <- permuted_results$distances
            return(permuted_distances[i, j])
          })
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        } else if (pvalue.method == "bootstrap") {
          # Bootstrap
          observed_distance <- distances[i, j]
          bootstrap_distances <- replicate(num.bootstraps, {
            # Extract a repetition sample
            bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
            bootstrap_data <- dataset[bootstrap_sample, ]
            bootstrap_results <- cmahalanobis(bootstrap_data, formula, plot = FALSE)
            bootstrap_distance <- bootstrap_results$distances[i, j]
            return(bootstrap_distance)
          })
          p_values[i, j] <- mean(abs(bootstrap_distances) >= abs(observed_distance))
        } else {
          stop("p_values method calculation not supported. Use 'chisq', 'permutation' or 'bootstrap'.")
        }
      }
    }
  }
  
  # Plot the heatmap if plot is TRUE
  if (plot) {
    requireNamespace("ggplot2")
    requireNamespace("reshape2")
    p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
    colnames(p_values_melt) <- c("Var1", "Var2", "value")
    
    print(ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "red") +
            ggplot2::labs(title = "p_values heatmap", x = "", y = "", fill = "p_value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  
  return(p_values)
}




#' Calculate Euclidean distance
#'
#' @name ceuclide
#' @title Calculate the Euclidean distance of a factor in a dataframe.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Euclidean distances about it.
#' @param dataset A dataframe.
#' @param formula The factor which you want to calculate the Euclidean distances matrix.
#' @param plot If TRUE, shows a plot of the Euclidean distances matrix.
#' @param plot_title The title of the plot.
#' @return The matrix containing distances.
#' @examples
#' 
#' # Example with iris dataset
#' 
#' ceuclide(iris, ~Species, plot = TRUE, plot_title = "Euclidean Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' ceuclide(mtcars, ~am, plot = TRUE, plot_title = "Euclidean Distance Between Groups")
#' 
#' @export
ceuclide <- function(dataset, formula, plot = TRUE, plot_title = "Euclidean Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variables in each group
  groups <- lapply(groups, function(df) {
    df <- df[, sapply(df, is.numeric)]  # Select only numeric columns
    return(df)
  })
  
  # Replace missing values with arithmetic mean in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means and covariances
  means <- lapply(groups, colMeans)
  covariances <- lapply(groups, function(df) {
    cov_matrix <- cov(df)
    n <- nrow(cov_matrix)
    reg <- 0.01  # Regularization value
    cov_matrix <- cov_matrix + diag(reg, nrow = n)
    return(cov_matrix)
  })
  
  # Create a empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Euclidean distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean of group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- sqrt(sum((mean_i - means[[j]])^2))
      }
    }
  }
  
  # If plot is TRUE, call the "plot_euclidean_distances" function and print the plot
  if (plot) {
    print(plot_euclidean_distances(distances, plot_title))
  }
  
  # Return a matrix containing distances 
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the mean 
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}


# Auxiliary function to print the Euclide distances plot 
plot_euclidean_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- reshape2::melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}


#' @name generate_report_ceuclide
#' @title Generate a Microsoft Word document about the Euclidean distance matrix and the p-values matrix with relative plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Euclidean distance matrix and the p-values matrix with relative plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Euclidean distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq".Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Euclidean distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ceuclide(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_ceuclide(mtcars, ~am)
#' 
#' @export
generate_report_ceuclide <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  ceuclide_results <- ceuclide(dataset, formula)
  distances <- ceuclide_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "chisq")  # Adjust method if needed
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "permutation")  # Adjust method if needed
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "bootstrap")  # Adjust method if needed
  }
  
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "reportceuclide.docx")
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Assuming the path to your template is correct
  rmarkdown::render(system.file("rmarkdown", "template_report_ceuclide.Rmd", package = "cmahalanobis"),
                    params = list(distances = distances, p_values = p_values),
                    output_file = "reportceuclide.docx")
}



#' @name pvaluesceucl
#' @title Calculate the p_values matrix for each species, using the Euclidean distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Euclidean distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Euclidean distances.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing the p_values matrix and, optionally, the plot.
#' #' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesceucl(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesceucl(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' 
#' @export
pvaluesceucl <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
  # Verify that input is a dataframe
  if (!is.data.frame(dataset)) {
    stop("dataset must be a dataframe")
  }
  
  # Extract the response variable name by the formula
  response <- all.vars(formula)[1]
  
  # Verify the presence of the response variable
  if (!response %in% names(dataset)) {
    stop(paste("The response variable", response, "isn't present"))
  }
  
  # Calculate Euclidean distances using "ceuclide" function
  ceuclide_results <- ceuclide(dataset, formula, plot = FALSE)
  distances <- ceuclide_results$distances
  
  # Obtain the groups number
  n <- nrow(distances)
  
  # Initialize p_value matrix
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate p_values
  for (i in 1:n) {
    df <- ncol(dataset) - 1  # Degrees of freedom (adjusted for matrix covariance estimate)
    for (j in 1:n) {
      if (i != j) {
        # Choose p_values method calculation based user inputs
        if (pvalue.method == "chisq") {
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE, log.p = TRUE)
        } else if (pvalue.method == "permutation") {
          # Permutation test
          observed_distance <- distances[i, j]
          permutation_distances <- replicate(num.permutations, {
            # Permute labels group and recalculate the distance
            permuted_data <- dataset
            permuted_data[[response]] <- sample(dataset[[response]])
            permuted_results <- ceuclide(permuted_data, formula, plot = FALSE)
            permuted_distances <- permuted_results$distances
            return(permuted_distances[i, j])
          })
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        } else if (pvalue.method == "bootstrap") {
          # Bootstrap
          observed_distance <- distances[i, j]
          bootstrap_distances <- replicate(num.bootstraps, {
            # Extract a sample with repetition 
            bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
            bootstrap_data <- dataset[bootstrap_sample, ]
            bootstrap_results <- ceuclide(bootstrap_data, formula, plot = FALSE)
            bootstrap_distance <- bootstrap_results$distances[i, j]
            return(bootstrap_distance)
          })
          p_values[i, j] <- mean(abs(bootstrap_distances) >= abs(observed_distance))
        } else {
          stop("p_values calculation method not supported. Use 'chisq', 'permutation' or 'bootstrap'.")
        }
      }
    }
  }
  
  # Create heatmap if plot is TRUE
  if (plot) {
    requireNamespace("reshape2")
    requireNamespace("ggplot2")
    p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
    colnames(p_values_melt) <- c("Var1", "Var2", "value")
    
    print(ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "red") +
            ggplot2::labs(title = "p_values heatmap", x = "", y = "", fill = "p_value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  
  return(p_values)
}



#' Calculate Manhattan distance
#' @name cmanhattan
#' @title Calculate a Manhattan distance of a factor in a dataframe.
#' @description This function takes a dataframe and a factor in input, and returns a matrix with the Manhattan distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Manhattan distance.
#' @param plot If TRUE, show a plot of distances.
#' @param plot_title The title of plot.
#' @return A matrix containing distances.
#' @examples
#' # Example with iris dataset
#' 
#' cmanhattan(iris, ~Species, plot = TRUE, plot_title = "Manhattan Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' cmanhattan(mtcars, ~am, plot = TRUE, plot_title = "Manhattan Distance Between Groups")
#' 
#' @export
cmanhattan <- function(dataset, formula, plot = TRUE, plot_title = "Manhattan Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variables in each group
  groups <- lapply(groups, function(df) {
    df <- df[, sapply(df, is.numeric)]  # Select only numeric columns
    return(df)
  })
  
  # Replace missing values with arithmetic mean in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means and covariances
  means <- lapply(groups, colMeans)
  covariances <- lapply(groups, function(df) {
    cov_matrix <- cov(df)
    n <- nrow(cov_matrix)
    reg <- 0.01  # Regularization value
    cov_matrix <- cov_matrix + diag(reg, nrow = n)
    return(cov_matrix)
  })
  
  # Create a empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Manhattan distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean in the group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- sum(abs(mean_i - means[[j]]))
      }
    }
  }
  
  # If plot is TRUE, call the "plot_manhattan_distances" function and print the plot
  if (plot) {
    print(plot_manhattan_distances(distances, plot_title))
  }
  
  # Return a matrix containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the mean 
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Imputazione multipla utilizzando mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Manhattan distances plot 
plot_manhattan_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- reshape2::melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}




#' @name generate_report_cmanhattan
#' @title Generate a Microsoft Word document about the Manhattan distance and the p-values matrices with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Manhattan distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Manhattan distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq".Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Manhattan distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cmanhattan(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cmanhattan(mtcars, ~am)
#' 
#' @export
generate_report_cmanhattan <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cmanhattan_results <- cmanhattan(dataset, formula)
  distances <- cmanhattan_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "chisq")  # Adjust method if needed
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "permutation")  # Adjust method if needed
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "bootstrap")  # Adjust method if needed
  }
  
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "reportcmanhattan.docx")
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Assuming the path to your template is correct
  rmarkdown::render(system.file("rmarkdown", "template_report_cmanhattan.Rmd", package = "cmahalanobis"),
                    params = list(distances = distances, p_values = p_values),
                    output_file = "reportcmanhattan.docx")
}

#' @name pvaluescmanh
#' @title Calculate the p_values matrix for each species, using Manhattan distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Manhattan distance for the distances calculation.
#' @param dataset A dataframe
#' @param formula A factor which you want to calculate Manhattan distances.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A matrix containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescmanh(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescmanh(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' 
#' @export
pvaluescmanh <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
  # Verify that input is a dataframe
  if (!is.data.frame(dataset)) {
    stop("dataset must be a dataframe")
  }
  
  # Extract the response variable name by the formula
  response <- all.vars(formula)[1]
  
  # Verify the presence of the response variable
  if (!response %in% names(dataset)) {
    stop(paste("The response variable", response, "isn't present"))
  }
  
  # Calculate Manhattan distances using "cmanhattan" function
  cmanhattan_results <- cmanhattan(dataset, formula, plot = FALSE)
  distances <- cmanhattan_results$distances
  
  # Obtain the groups number
  n <- nrow(distances)
  
  # Initialize p_value matrix
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate p_values
  for (i in 1:n) {
    df <- ncol(dataset) - 1  # Degrees of freedom (adjusted for matrix covariance estimate)
    for (j in 1:n) {
      if (i != j) {
        # Choose p_values method calculation based user inputs
        if (pvalue.method == "chisq") {
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE, log.p = TRUE)
        } else if (pvalue.method == "permutation") {
          # Permutation test
          observed_distance <- distances[i, j]
          permutation_distances <- replicate(num.permutations, {
            # Permute labels group and recalculate the distance
            permuted_data <- dataset
            permuted_data[[response]] <- sample(dataset[[response]])
            permuted_results <- cmanhattan(permuted_data, formula, plot = FALSE)
            permuted_distances <- permuted_results$distances
            return(permuted_distances[i, j])
          })
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        } else if (pvalue.method == "bootstrap") {
          # Bootstrap
          observed_distance <- distances[i, j]
          bootstrap_distances <- replicate(num.bootstraps, {
            # Extract a sample with repetition 
            bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
            bootstrap_data <- dataset[bootstrap_sample, ]
            bootstrap_results <- cmanhattan(bootstrap_data, formula, plot = FALSE)
            bootstrap_distance <- bootstrap_results$distances[i, j]
            return(bootstrap_distance)
          })
          p_values[i, j] <- mean(abs(bootstrap_distances) >= abs(observed_distance))
        } else {
          stop("p_values calculation method not supported. Use 'chisq', 'permutation' or 'bootstrap'.")
        }
      }
    }
  }
  
  # Create heatmap if plot is TRUE
  if (plot) {
    requireNamespace("reshape2")
    requireNamespace("ggplot2")
    p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
    colnames(p_values_melt) <- c("Var1", "Var2", "value")
    
    print(ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "red") +
            ggplot2::labs(title = "p_values heatmap", x = "", y = "", fill = "p_value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  
  return(p_values)
}

#' @name cchebyshev
#' @title Calculate the p_values matrix for each species, using Chebyshev distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using the Chebyshev distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Chebyshev distance.
#' @param plot If TRUE, displays a plot of distances.
#' @param plot_title The title of plot.
#' @return A matrix containing distances and, optionally, the plot.
#' @examples 
#' # Example with iris dataset
#' 
#' cchebyshev(iris, ~Species, plot = TRUE, plot_title = "Chebyshev Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' cchebyshev(mtcars, ~am, plot = TRUE, plot_title = "Chebyshev Distance Between Groups")
#' 
#' 
#' @export
cchebyshev <- function(dataset, formula, plot = TRUE, plot_title = "Chebyshev Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variables in each group
  groups <- lapply(groups, function(df) {
    df <- df[, sapply(df, is.numeric)]  # Select only numeric columns
    return(df)
  })
  
  # Replace missing values with arithmetic mean in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means and covariances
  means <- lapply(groups, colMeans)
  covariances <- lapply(groups, function(df) {
    cov_matrix <- cov(df)
    n <- nrow(cov_matrix)
    reg <- 0.01  # Regularization value
    cov_matrix <- cov_matrix + diag(reg, nrow = n)
    return(cov_matrix)
  })
  
  # Create a empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Chebyshev distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]]  # Precalculated mean in the group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- max(abs(mean_i - means[[j]]))
      }
    }
  }
  
  # If plot is TRUE, call the "plot_chebyshev_distances" function and print the plot
  if (plot) {
    print(plot_chebyshev_distances(distances, plot_title))
  }
  
  # Return a list containing distances 
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the mean 
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}


# Auxiliary function to print the Mahalanobis distances plot 
plot_chebyshev_distances <- function(distances, plot_title) {
  requireNamespace("ggplot2")
  requireNamespace("reshape2")
  output_df <- as.data.frame(distances)
  output_df$Species <- rownames(output_df)
  
  output_df_long <- reshape2::melt(output_df, id.vars = "Species")
  colnames(output_df_long) <- c("Species", "Comparison", "Distance")
  
  ggplot2::ggplot(output_df_long, ggplot2::aes(x = Species, y = Distance, fill = Comparison)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(title = plot_title, x = "Species", y = "Distance", fill = "Comparison") +
    ggplot2::theme_minimal()
}


#' @name generate_report_cchebyshev
#' @title Generate a Microsoft Word document about the Chebyshev distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Chebyshev distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Chebyshev distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Chebyshev distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cchebyshev(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cchebyshev(mtcars, ~am)
#' 
#' @export
generate_report_cchebyshev <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cchebyshev_results <- cchebyshev(dataset, formula)
  distances <- cchebyshev_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "chisq")  # Adjust method if needed
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "permutation")  # Adjust method if needed
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "bootstrap")  # Adjust method if needed
  }
  
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "reportcchebyshev.docx")
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  
  
  # Assuming the path to your template is correct
  rmarkdown::render(system.file("rmarkdown", "template_report_cchebyshev.Rmd", package = "cmahalanobis"),
                    params = list(distances = distances, p_values = p_values),
                    output_file = "reportcchebyshev.docx")
}





#' @name pvaluesccheb
#' @title Calculate the p_values matrix for each species, using Chebyshev distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Chebyshev distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Chebyshev distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccheb(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccheb(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluesccheb <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
  # Verify that input is a dataframe
  if (!is.data.frame(dataset)) {
    stop("dataset must be a dataframe")
  }
  
  # Extract the response variable name by the formula
  response <- all.vars(formula)[1]
  
  # Verify the presence of the response variable
  if (!response %in% names(dataset)) {
    stop(paste("The response variable", response, "isn't present"))
  }
  
  # Calculate Chebyshev distances using "cchebyshev" function
  cchebyshev_results <- cchebyshev(dataset, formula, plot = FALSE)
  distances <- cchebyshev_results$distances
  
  # Obtain the groups number
  n <- nrow(distances)
  
  # Initialize p_value matrix
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate p_values
  for (i in 1:n) {
    df <- ncol(dataset) - 1  # Degrees of freedom (adjusted for matrix covariance estimate)
    for (j in 1:n) {
      if (i != j) {
        # Choose p_values method calculation based user inputs
        if (pvalue.method == "chisq") {
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE, log.p = TRUE)
        } else if (pvalue.method == "permutation") {
          # Permutation test
          observed_distance <- distances[i, j]
          permutation_distances <- replicate(num.permutations, {
            # Permute labels group and recalculate the distance
            permuted_data <- dataset
            permuted_data[[response]] <- sample(dataset[[response]])
            permuted_results <- cchebyshev(permuted_data, formula, plot = FALSE)
            permuted_distances <- permuted_results$distances
            return(permuted_distances[i, j])
          })
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        } else if (pvalue.method == "bootstrap") {
          # Bootstrap
          observed_distance <- distances[i, j]
          bootstrap_distances <- replicate(num.bootstraps, {
            # Extract a sample with repetition 
            bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
            bootstrap_data <- dataset[bootstrap_sample, ]
            bootstrap_results <- cchebyshev(bootstrap_data, formula, plot = FALSE)
            bootstrap_distance <- bootstrap_results$distances[i, j]
            return(bootstrap_distance)
          })
          p_values[i, j] <- mean(abs(bootstrap_distances) >= abs(observed_distance))
        } else {
          stop("p_values calculation method not supported. Use 'chisq', 'permutation' or 'bootstrap'.")
        }
      }
    }
  }
  
  # Create heatmap if plot is TRUE
  if (plot) {
    requireNamespace("reshape2")
    requireNamespace("ggplot2")
    p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
    colnames(p_values_melt) <- c("Var1", "Var2", "value")
    
    print(ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "red") +
            ggplot2::labs(title = "p_values heatmap", x = "", y = "", fill = "p_value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  
  return(p_values)
}

