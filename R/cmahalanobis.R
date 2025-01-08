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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
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
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Mahalanobis distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean of the group i
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
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescmaha(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcmahalanobis_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cmahalanobis.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcmahalanobis.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcmahalanobis_files", recursive = TRUE)
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
  # Verify that input is a dataframe
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response variable name by the formula
  response <- all.vars(formula)[1]
  
  # Verify the presence of the response variable
  if (!response %in% names(dataset)) {
    stop(paste("Response variable", response, "is absent in the dataset"))
  }
  
  # Calculate Mahalanobis distances
  mahalanobis_results <- cmahalanobis(dataset, formula, plot = FALSE)
  distances <- mahalanobis_results$distances
  
  # Obtain the number of groups
  n <- nrow(distances)
  
  # Initialize p_value matrix
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate p_values
  for (i in 1:n) {
    df <- ncol(dataset) - 1  # Degrees of freedom (adjusted for covariance matrix estimate)
    for (j in 1:n) {
      if (i != j) {
        if (pvalue.method == "chisq") {
          p_values[i, j] <- pchisq(distances[i, j], df, lower.tail = FALSE, log.p = TRUE)
        } else if (pvalue.method == "permutation") {
          observed_distance <- distances[i, j]
          permutation_distances <- replicate(num.permutations, {
            permuted_data <- dataset
            permuted_data[[response]] <- sample(dataset[[response]])
            permuted_results <- cmahalanobis(permuted_data, formula, plot = FALSE)
            permuted_distances <- permuted_results$distances
            return(permuted_distances[i, j])
          })
          p_values[i, j] <- mean(permutation_distances >= observed_distance)
        } else if (pvalue.method == "bootstrap") {
          observed_distance <- distances[i, j]
          bootstrap_distances <- replicate(num.bootstraps, {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means 
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Euclidean distance between each pair of groups
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

# Auxiliary function to impute missing values with the multiple imputation 
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
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesceucl(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportceuclide_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_ceuclide.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportceuclide.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportceuclide_files", recursive = TRUE)
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means 
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Manhattan distance between each pair of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean of group i
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

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
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
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescmanh(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcmanhattan_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cmanhattan.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcmanhattan.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcmanhattan_files", recursive = TRUE)
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Chebyshev distance between each pair of groups
  for (i in 1:n) {
    mean_i <- means[[i]]  # Precalculated mean of group i
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

# Auxiliary function to impute missing values with the multiple imputation
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
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesccheb(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcchebyshev_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cchebyshev.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcchebyshev.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcchebyshev_files", recursive = TRUE)
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



#' Calculate Hamming distance
#'
#' @name chamming
#' @title Calculate the Hamming distance of a factor in a dataframe.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Hamming distances about it.
#' @param dataset A dataframe.
#' @param formula The factor which you want to calculate the Hamming distances matrix.
#' @param plot If TRUE, shows a plot of the Hamming distances matrix.
#' @param plot_title The title of the plot.
#' @return The matrix containing distances.
#' @examples
#' 
#' # Example with iris dataset
#' 
#' chamming(iris, ~Species, plot = TRUE, plot_title = "Hamming Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' chamming(mtcars, ~am, plot = TRUE, plot_title = "Hamming Distance Between Groups")
#' 
#' @export
chamming <- function(dataset, formula, plot = TRUE, plot_title = "Hamming Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variable in each group
  groups <- lapply(groups, function(df) {
    df <- df[, sapply(df, is.numeric)] # Select only numeric columns
    return(df)
  })
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Hamming distance between each couple of groups
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- mean(apply(groups[[i]], 1, function(row_i) {
          apply(groups[[j]], 1, function(row_j) {
            sum(row_i != row_j)
          })
        }))
      }
    }
  }
  
  # If plot is TRUE, call the "plot_hamming_distances" function and print the plot
  if (plot) {
    print(plot_hamming_distances(distances, plot_title))
  }
  
  # Return a list containing distances 
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Hamming distances plot
plot_hamming_distances <- function(distances, plot_title) {
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


#' @name pvalueschamm
#' @title Calculate the p_values matrix for each species, using Hamming distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Hamming distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Hamming distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvalueschamm(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvalueschamm(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvalueschamm <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Hamming distances using "chamming" function
  chamming_results <- chamming(dataset, formula, plot = FALSE)
  distances <- chamming_results$distances
  
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
            permuted_results <- chamming(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- chamming(bootstrap_data, formula, plot = FALSE)
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

#' @name generate_report_chamming
#' @title Generate a Microsoft Word document about the Hamming distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Hamming distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Hamming distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Hamming distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_chamming(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_chamming(mtcars, ~am)
#' 
#' @export
generate_report_chamming <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  chamming_results <- chamming(dataset, formula)
  distances <- chamming_results
  if (pvalue.method == "chisq") {
    p_values <- pvalueschamm(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvalueschamm(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvalueschamm(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportchamming_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_chamming.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportchamming.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportchamming_files", recursive = TRUE)
}


#' @name ccanberra
#' @title Calculate the Canberra distance for each species.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Canberra distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Canberra distances matrix.
#' @param plot Logical, if TRUE, a plot of Canberra distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Canberra distances between each pair of groups and the plot.
#'
#'
#' @examples
#' # Example with the iris dataset
#'
#' ccanberra(iris, ~Species, plot = TRUE, plot_title = "Canberra Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' ccanberra(mtcars, ~am, plot = TRUE, plot_title = "Canberra Distance Between Groups")
#' 
#' @export
ccanberra <- function(dataset, formula, plot = TRUE, plot_title = "Canberra Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate the Canberra distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean of group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- mean(apply(groups[[j]], 1, function(row) sum(abs(row - mean_i) / (abs(row) + abs(mean_i)))))
      }
    }
  }
  
  # If plot is TRUE, call the "plot_canberra_distances" and print the plot
  if (plot) {
    print(plot_canberra_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}


# Auxiliary function to impute missing values with the multiple imputation 
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Canberra distances plot
plot_canberra_distances <- function(distances, plot_title) {
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



#' @name pvaluesccanb
#' @title Calculate the p_values matrix for each species, using Canberra distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Canberra distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Canberra distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccanb(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccanb(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluesccanb <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Canberra distances using "ccanberra" function
  ccanberra_results <- ccanberra(dataset, formula, plot = FALSE)
  distances <- ccanberra_results$distances
  
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
            permuted_results <- ccanberra(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- ccanberra(bootstrap_data, formula, plot = FALSE)
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


#' @name generate_report_ccanberra
#' @title Generate a Microsoft Word document about the Canberra distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Canberra distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Canberra distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Canberra distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ccanberra(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_ccanberra(mtcars, ~am)
#' 
#' @export
generate_report_ccanberra <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  ccanberra_results <- ccanberra(dataset, formula)
  distances <- ccanberra_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluesccanb(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesccanb(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesccanb(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportccanberra_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_ccanberra.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportccanberra.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportccanberra_files", recursive = TRUE)
}


#' Calculate Minkowski distance
#'
#' @name cminkowski
#' @title Calculate the Minkowski distance of a factor in a dataframe.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Minkowski distances about it.
#' @param dataset A dataframe.
#' @param formula The factor which you want to calculate the Minkowski distances matrix.
#' @param p Order of the Minkowski distance
#' @param plot If TRUE, shows a plot of the Minkowski distances matrix.
#' @param plot_title The title of the plot.
#' @return The matrix containing distances.
#' @examples
#' 
#' # Example with iris dataset
#' 
#' cminkowski(iris, ~Species, p = 3, plot = TRUE, plot_title = "Minkowski Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' cminkowski(mtcars, ~am, p = 3, plot = TRUE, plot_title = "Minkowski Distance Between Groups")
#' 
#' @export
cminkowski <- function(dataset, formula, p = 3, plot = TRUE, plot_title = "Minkowski Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Minkowski distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean of group i
    for (j in 1:n) {
      if (i != j) {
        distances[i, j] <- mean(apply(groups[[j]], 1, function(row) {
          sum(abs(row - mean_i)^p)^(1/p)
        }))
      }
    }
  }
  
  # If plot is TRUE, call the function "plot_minkowski_distances" and print the plot
  if (plot) {
    print(plot_minkowski_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Minkowski distances plot
plot_minkowski_distances <- function(distances, plot_title) {
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


#' @name pvaluescmink
#' @title Calculate the p_values matrix for each species, using Minkowski distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Minkowski distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Minkowski distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param p Order of the Minkowski distance
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescmink(iris,~Species, p = 3, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescmink(mtcars,~am, p = 3, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescmink <- function(dataset, formula, p = 3, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Minkowski distances using "cminkowski" function
  cminkowski_results <- cminkowski(dataset, formula, p = p, plot = FALSE)
  distances <- cminkowski_results$distances
  
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
            permuted_results <- cminkowski(permuted_data, formula, p = p, plot = FALSE)
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
            bootstrap_results <- cminkowski(bootstrap_data, formula, p = p, plot = FALSE)
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

#' @name generate_report_cminkowski
#' @title Generate a Microsoft Word document about the Minkowski distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Minkowski distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Minkowski distance matrix and the p_values matrix.
#' @param p Order of the Minkowski distance
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Minkowski distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cminkowski(iris, ~Species, p = 3)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cminkowski(mtcars, ~am, p = 3)
#' 
#' @export
generate_report_cminkowski <- function(dataset, formula, p = 3, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cminkowski_results <- cminkowski(dataset, formula, p = p)
  distances <- cminkowski_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescmink(dataset, formula, p = p, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescmink(dataset, formula, p = p, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescmink(dataset, formula, p = p, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcminkowski_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cminkowski.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcminkowski.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcminkowski_files", recursive = TRUE)
}


#' Calculate Cosine distance
#'
#' @name ccosine
#' @title Calculate the Cosine distance of a factor in a dataframe.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Cosine distances about it.
#' @param dataset A dataframe.
#' @param formula The factor which you want to calculate the Cosine distances matrix.
#' @param plot If TRUE, shows a plot of the Cosine distances matrix.
#' @param plot_title The title of the plot.
#' @return The matrix containing distances.
#' @examples
#' 
#' # Example with iris dataset
#' 
#' ccosine(iris, ~Species, plot = TRUE, plot_title = "Cosine Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' ccosine(mtcars, ~am, plot = TRUE, plot_title = "Cosine Distance Between Groups")
#' 
#' @export
ccosine <- function(dataset, formula, plot = TRUE, plot_title = "Cosine Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Cosine distance between each couple of groups
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        num <- sum(means[[i]] * means[[j]])
        denom <- sqrt(sum(means[[i]]^2)) * sqrt(sum(means[[j]]^2))
        distances[i, j] <- 1 - (num / denom)
      }
    }
  }
  
  
  # If plot is TRUE, call the "plot_cosine_distances" function and print the plot
  if (plot) {
    print(plot_cosine_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Cosine distances plot
plot_cosine_distances <- function(distances, plot_title) {
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



#' @name pvaluesccosi
#' @title Calculate the p_values matrix for each species, using Cosine distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Cosine distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Cosine distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccosi(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccosi(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluesccosi <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Cosine distances using "ccosine" function
  ccosine_results <- ccosine(dataset, formula, plot = FALSE)
  distances <- ccosine_results$distances
  
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
            permuted_results <- ccosine(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- ccosine(bootstrap_data, formula, plot = FALSE)
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



#' @name generate_report_ccosine
#' @title Generate a Microsoft Word document about the Cosine distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Cosine distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Cosine distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Cosine distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ccosine(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_ccosine(mtcars, ~am)
#' 
#' @export
generate_report_ccosine <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  ccosine_results <- ccosine(dataset, formula)
  distances <- ccosine_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluesccosi(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluesccosi(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluesccosi(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportccosine_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_ccosine.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportccosine.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportccosine_files", recursive = TRUE)
}


#' @name ccbhattacharyya
#' @title Calculate the Bhattacharyya distance for each species.
#'
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Bhattacharyya distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Bhattacharyya distances matrix.
#' @param plot Logical, if TRUE, a plot of Bhattacharyya distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Bhattacharyya distances between each pair of groups and the plot.
#' @examples
#' # Example with the iris dataset
#' cbhattacharyya(iris, ~Species, plot = TRUE, plot_title = "Bhattacharyya Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' cbhattacharyya(mtcars, ~am, plot = TRUE, plot_title = "Bhattacharyya Distance Between Groups")
#' 
#' @export
cbhattacharyya <- function(dataset, formula, plot = TRUE, plot_title = "Bhattacharyya Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
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
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Bhattacharyya distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean in the group i
    cov_i <- covariances[[i]]
    for (j in 1:n) {
      if (i != j) {
        mean_j <- means[[j]]
        cov_j <- covariances[[j]]
        
        # Calculate the average covariance matrix
        cov_mean <- (cov_i + cov_j) / 2
        
        # Calculate the Bhattacharyya distance
        term1 <- 0.25 * t(mean_i - mean_j) %*% solve(cov_mean) %*% (mean_i - mean_j)
        term2 <- 0.5 * log(det(cov_mean) / sqrt(det(cov_i) * det(cov_j)))
        distances[i, j] <- term1 + term2
      }
    }
  }
  
  # If plot is TRUE, call the "plot_bhattacharyya_distances" function and print the plot
  if (plot) {
    print(plot_bhattacharyya_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Bhattacharyya distances plot
plot_bhattacharyya_distances <- function(distances, plot_title) {
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

#' @name pvaluescbatt
#' @title Calculate the p_values matrix for each species, using Bhattacharyya distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Bhattacharyya distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Bhattacharyya distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescbatt(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescbatt(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescbatt <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Bhattacharyya distances using "cbhattacharyya" function
  cbhattacharyya_results <- cbhattacharyya(dataset, formula, plot = FALSE)
  distances <- cbhattacharyya_results$distances
  
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
            permuted_results <- cbhattacharyya(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- cbhattacharyya(bootstrap_data, formula, plot = FALSE)
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


#' @name generate_report_cbhattacharyya
#' @title Generate a Microsoft Word document about the Bhattacharyya distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Bhattacharyya distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Bhattacharyya distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Bhattacharyya distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cbhattacharyya(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cbhattacharyya(mtcars, ~am)
#' 
#' @export
generate_report_cbhattacharyya <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cbhattacharyya_results <- cbhattacharyya(dataset, formula)
  distances <- cbhattacharyya_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescbatt(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescbatt(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescbatt(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcbhattacharyya_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cbhattacharyya.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcbhattacharyya.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcbhattacharyya_files", recursive = TRUE)
}


#' @name cjaccard
#' @title Calculate the Jaccard distance for each species.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Jaccard distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Jaccard distances matrix.
#' @param plot Logical, if TRUE, a plot of Jaccard distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Jaccard distances between each pair of groups and the plot.
#'
#'
#' @examples
#' # Example with the iris dataset
#'
#' cjaccard(iris, ~Species, plot = TRUE, plot_title = "Jaccard Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' cjaccard(mtcars, ~am, plot = TRUE, plot_title = "Jaccard Distance Between Groups")
#' 
#' @export
cjaccard <- function(dataset, formula, plot = TRUE, plot_title = "Jaccard Distance Between Groups") {
  # Verify that the input is a data frame
  if (!is.data.frame(dataset)) {
    stop("The input must be a dataframe")
  }
  
  # Extract the response and predictor variables from the formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  
  # Split the data into groups based on the response variable
  groups <- split(dataset, dataset[[response]])
  
  # Select only numeric variables in each group and binarize them
  binarize <- function(df) {
    df <- df[, sapply(df, is.numeric)]  # Select only numeric columns
    df <- as.data.frame(lapply(df, function(x) as.integer(x > mean(x))))
    return(df)
  }
  
  groups <- lapply(groups, binarize)
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Jaccard distance between each couple of groups
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Flatten groups to create binary vectors
        vector_i <- as.vector(unlist(groups[[i]]))
        vector_j <- as.vector(unlist(groups[[j]]))
        
        # Calculate Jaccard similarity
        intersection <- sum(vector_i & vector_j)
        union <- sum(vector_i | vector_j)
        jaccard_similarity <- intersection / union
        
        # Calculate Jaccard distance
        distances[i, j] <- 1 - jaccard_similarity
      }
    }
  }
  
  # If plot is TRUE, call the "plot_jaccard_distances" function and print the plot
  if (plot) {
    print(plot_jaccard_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Jaccard distances plot
plot_jaccard_distances <- function(distances, plot_title) {
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



#' @name pvaluescjacc
#' @title Calculate the p_values matrix for each species, using Jaccard distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Jaccard distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Jaccard distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescjacc(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescjacc(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescjacc <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Jaccard distances using "cjaccard" function
  cjaccard_results <- cjaccard(dataset, formula, plot = FALSE)
  distances <- cjaccard_results$distances
  
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
            permuted_results <- cjaccard(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- cjaccard(bootstrap_data, formula, plot = FALSE)
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


#' @name generate_report_cjaccard
#' @title Generate a Microsoft Word document about the Jaccard distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Jaccard distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Jaccard distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Jaccard distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cjaccard(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cjaccard(mtcars, ~am)
#' 
#' @export
generate_report_cjaccard <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cjaccard_results <- cjaccard(dataset, formula)
  distances <- cjaccard_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescjacc(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescjacc(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescjacc(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcjaccard_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cjaccard.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcjaccard.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcjaccard_files", recursive = TRUE)
}


#' @name chellinger
#' @title Calculate the Hellinger distance for each species.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Hellinger distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Hellinger distances matrix.
#' @param plot Logical, if TRUE, a plot of Hellinger distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Hellinger distances between each pair of groups and the plot.
#' @examples
#' # Example with the iris dataset
#'
#' chellinger(iris, ~Species, plot = TRUE, plot_title = "Hellinger Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' chellinger(mtcars, ~am, plot = TRUE, plot_title = "Hellinger Distance Between Groups")
#' 
#' @export
chellinger <- function(dataset, formula, plot = TRUE, plot_title = "Hellinger Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Precalculate means and standard deviations
  means <- lapply(groups, colMeans)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Hellinger distance between each couple of groups
  for (i in 1:n) {
    mean_i <- means[[i]] # Precalculated mean in group i
    for (j in 1:n) {
      if (i != j) {
        mean_j <- means[[j]]
        
        # Calculate Hellinger distance
        hellinger_distance <- sqrt(0.5 * sum((sqrt(mean_i) - sqrt(mean_j))^2))
        distances[i, j] <- hellinger_distance
      }
    }
  }
  
  # If plot is TRUE, call the "plot_hellinger_distances" function and print the plot
  if (plot) {
    print(plot_hellinger_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
  return(imputedData)
}

# Auxiliary function to print the Hellinger distances plot
plot_hellinger_distances <- function(distances, plot_title) {
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

#' @name pvalueschell
#' @title Calculate the p_values matrix for each species, using Hellinger distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Hellinger distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Hellinger distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvalueschell(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvalueschell(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvalueschell <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Hellinger distances using "chellinger" function
  chellinger_results <- chellinger(dataset, formula, plot = FALSE)
  distances <- chellinger_results$distances
  
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
            permuted_results <- chellinger(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- chellinger(bootstrap_data, formula, plot = FALSE)
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




#' @name generate_report_chellinger
#' @title Generate a Microsoft Word document about the Hellinger distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Hellinger distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Hellinger distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Hellinger distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_chellinger(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_chellinger(mtcars, ~am)
#' 
#' @export
generate_report_chellinger <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  chellinger_results <- chellinger(dataset, formula)
  distances <- chellinger_results
  if (pvalue.method == "chisq") {
    p_values <- pvalueschell(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvalueschell(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvalueschell(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportchellinger_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_chellinger.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportchellinger.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportchellinger_files", recursive = TRUE)
}



#' @name cbraycurtis
#' @title Calculate the Bray-Curtis distance for each species.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Bray-Curtis distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Bray-Curtis distances matrix.
#' @param plot Logical, if TRUE, a plot of Bray-Curtis distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Bray-Curtis distances between each pair of groups and the plot.
#' @examples
#' # Example with the iris dataset
#'
#' cbraycurtis(iris, ~Species, plot = TRUE, plot_title = "Bray-Curtis Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' cbraycurtis(mtcars, ~am, plot = TRUE, plot_title = "Bray-Curtis Distance Between Groups")
#' 
#' @export
cbraycurtis <- function(dataset, formula, plot = TRUE, plot_title = "Bray-Curtis Distance Between Groups") {
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
  })
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Bray-Curtis distance between each couple of groups
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        common_columns <- intersect(colnames(groups[[i]]), colnames(groups[[j]]))
        group_i <- rowSums(groups[[i]][, common_columns, drop = FALSE])
        group_j <- rowSums(groups[[j]][, common_columns, drop = FALSE])
        
        sum_abs_min <- sum(abs(group_i - group_j))
        
        sum_total <- sum(group_i) + sum(group_j)
        
        bray_curtis_distance <- sum_abs_min / sum_total
        distances[i, j] <- bray_curtis_distance
      }
    }
  }
  
  
  # If plot is TRUE, call the "plot_braycurtis_distances" function and print the plot
  if (plot) {
    print(plot_braycurtis_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
}

# Auxiliary function to print the Bray-Curtis distances plot
plot_braycurtis_distances <- function(distances, plot_title) {
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





#' @name pvaluescbrcu
#' @title Calculate the p_values matrix for each species, using Bray-Curtis distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Bray-Curtis distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Bray-Curtis distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescbrcu(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescbrcu(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescbrcu <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Bray-Curtis distances using "cbraycurtis" function
  cbraycurtis_results <- cbraycurtis(dataset, formula, plot = FALSE)
  distances <- cbraycurtis_results$distances
  
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
            permuted_results <- cbraycurtis(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- cbraycurtis(bootstrap_data, formula, plot = FALSE)
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

#' @name generate_report_cbraycurtis
#' @title Generate a Microsoft Word document about the Bray-Curtis distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Bray-Curtis distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Bray-Curtis distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Bray-Curtis distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cbraycurtis(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cbraycurtis(mtcars, ~am)
#' 
#' @export
generate_report_cbraycurtis <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  cbraycurtis_results <- cbraycurtis(dataset, formula)
  distances <- cbraycurtis_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescbrcu(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescbrcu(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescbrcu(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcbraycurtis_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_cbraycurtis.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcbraycurtis.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcbraycurtis_files", recursive = TRUE)
}


#' @name csorensendice
#' @title Calculate the Sorensen-Dice distance for each species.
#' @description
#' This function takes a dataframe and a factor in input, and returns a matrix with the Sorensen-Dice distances about it.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Sorensen-Dice distances matrix.
#' @param plot Logical, if TRUE, a plot of Sorensen-Dice distances matrix is displayed.
#' @param plot_title The title to be used for the plot if plot is TRUE.
#' @return A matrix containing Sorensen-Dice distances between each pair of groups and the plot.
#' @examples
#' # Example with the iris dataset
#'
#' csorensendice(iris, ~Species, plot = TRUE, plot_title = "Sorensen-Dice Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' csorensendice(mtcars, ~am, plot = TRUE, plot_title = "Sorensen-Dice Distance Between Groups")
#' 
#' @export
csorensendice <- function(dataset, formula, plot = TRUE, plot_title = "Sorensen-Dice Distance Between Groups") {
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
  
  # Replace missing values with the multiple imputation in each dataframe into the list
  groups <- lapply(groups, impute_with_multiple_imputation)
  
  # Obtain the number of groups
  n <- length(groups)
  
  # Create an empty matrix to store distances
  distances <- matrix(0, nrow = n, ncol = n)
  
  # Calculate Sorensen-Dice distance between each couple of groups
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        common_columns <- intersect(colnames(groups[[i]]), colnames(groups[[j]]))
        group_i <- rowSums(groups[[i]][, common_columns, drop = FALSE])
        group_j <- rowSums(groups[[j]][, common_columns, drop = FALSE])
        
        intersection_sum <- sum(pmin(group_i, group_j))
        
        sorensen_dice_distance <- (2 * intersection_sum) / (sum(group_i) + sum(group_j))
        distances[i, j] <- sorensen_dice_distance
      }
    }
  }
  
  # If plot is TRUE, call the "plot_sorensendice_distances" function and print the plot
  if (plot) {
    print(plot_sorensendice_distances(distances, plot_title))
  }
  
  # Return a list containing distances
  return(list(distances = distances))
}

# Auxiliary function to impute missing values with the multiple imputation
impute_with_multiple_imputation <- function(df, m = 5, method = "cart") {
  # Multiple imputation using mice
  imp <- mice(df, m = m, method = method)
  imputedData <- complete(imp, action = 1)
}

# Auxiliary function to print the Sorensen-Dice distances plot
plot_sorensendice_distances <- function(distances, plot_title) {
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



#' @name pvaluescsore
#' @title Calculate the p_values matrix for each species, using Sorensen-Dice distance as a base.
#' @description
#' This function takes a dataset, a factor, a p_value method, number of bootstraps and permutation when necessary, and returns a p_values matrix between each pair of species and a plot if the user select TRUE using Sorensen-Dice distance for the distances calculation.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate Sorensen-Dice distance.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @param plot if TRUE, plot the p_values heatmap. The default value is TRUE.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescsore(iris,~Species, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescsore(mtcars,~am, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10)
#' @export
pvaluescsore <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 100, num.bootstraps = 10, plot = TRUE) {
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
  
  # Calculate Sorensen-Dice distances using "csorensendice" function
  csorensendice_results <- csorensendice(dataset, formula, plot = FALSE)
  distances <- csorensendice_results$distances
  
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
            permuted_results <- csorensendice(permuted_data, formula, plot = FALSE)
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
            bootstrap_results <- csorensendice(bootstrap_data, formula, plot = FALSE)
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

#' @name generate_report_csorensendice
#' @title Generate a Microsoft Word document about the Sorensen-Dice distance matrix and the p-values matrix with corresponding plots.
#'
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Sorensen-Dice distance matrix and the p-values matrix with corresponding plots.
#' @param dataset A dataframe.
#' @param formula A factor which you want to calculate the Sorensen-Dice distance matrix and the p_values matrix.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "chisq". Other methods are "permutation" and "bootstrap".
#' @param num.permutations Number of permutation to specify if you select "permutation" in "pvalue.method". The default value is 100.
#' @param num.bootstraps Number of bootstrap to specify if you select "bootstrap" in "p_value method". The default value is 10.
#' @return A Microsoft Word document about the Sorensen-Dice distance matrix and the p_values matrix.
#' @examples
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_csorensendice(iris, ~Species)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_csorensendice(mtcars, ~am)
#' 
#' @export
generate_report_csorensendice <- function(dataset, formula, pvalue.method = "chisq", num.permutations = 10, num.bootstraps = 10) {
  requireNamespace("rmarkdown")
  csorensendice_results <- csorensendice(dataset, formula)
  distances <- csorensendice_results
  if (pvalue.method == "chisq") {
    p_values <- pvaluescsore(dataset, formula, pvalue.method = "chisq")
  } else if (pvalue.method == "permutation") {
    p_values <- pvaluescsore(dataset, formula, pvalue.method = "permutation")
  } else if (pvalue.method == "bootstrap") {
    p_values <- pvaluescsore(dataset, formula, pvalue.method = "bootstrap")
  }
  
  dir_path <- file.path(getwd(), "reportcsorensendice_files/figure-docx")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  template_path <- system.file("rmarkdown", "template_report_csorensendice.Rmd", package = "cmahalanobis")
  if (file.exists(template_path)) {
    message("Template found at: ", template_path)
  } else {
    stop("Template not found!")
  }
  
  output_file <- "reportcsorensendice.docx"
  tryCatch({
    rmarkdown::render(template_path,
                      params = list(distances = distances, p_values = p_values),
                      output_file = output_file)
  }, error = function(e) {
    message("Error during rendering: ", e$message)
  })
  
  get_desktop_path <- function() {
    home <- Sys.getenv("HOME")
    sysname <- Sys.info()["sysname"]
    if (sysname == "Windows") {
      return(file.path(Sys.getenv("USERPROFILE"), "Desktop"))
    } else if (sysname == "Darwin") {  
      return(file.path(home, "Desktop"))
    } else if (sysname == "Linux") {  
      return(file.path(home, "Desktop"))
    } else {
      stop("Unsupported OS")
    }
  }
  
  desktop_path <- get_desktop_path()
  file.rename(output_file, file.path(desktop_path, output_file))
  
  unlink("reportcsorensendice_files", recursive = TRUE)
}
