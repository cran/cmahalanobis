# cmahalanobis.R

utils::globalVariables(c("Var1", "Var2", "value"))

#' @importFrom stats mahalanobis
#' @importFrom stats cov
#' @importFrom stats model.matrix
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 aes
#' @importFrom reshape2 melt 
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats dist
#' @importFrom stats na.omit
#' @importFrom stats pchisq
#' @importFrom stats pchisq
#' @importFrom stats sd
#' @importFrom gridExtra grid.arrange
#' @importFrom matrixStats colSums2 colMeans2 sum2 mean2 rowSums2
#' 
#' @name cmahalanobis
#' @title Calculate the Mahalanobis distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Mahalanobis distances about each pair of factors inside them. You can also select "index" to calculate the Mahalanobis distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Mahalanobis distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Mahalanobis distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore factors, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @param pvalues_chisq If TRUE, print the result of the chi-squared test on squared distances. The distances with "pvalues_chisq = FALSE" are not squared; instead, with "pvalues_chisq = TRUE", the squared Mahalanobis distances with corresponding p_values will be printed. Default is FALSE.
#' @return According to the option chosen in formula and in pvalues_chisq, with "index" and "pvalues_chisq = TRUE" the squared Mahalanobis distance matrix will be printed with corresponding pvalues; instead, with "index" and "pvalues_chisq = FALSE", only the Mahalanobis distances (not squared) will be printed. By specifying variables, the Mahalanobis distances matrix or matrices (two or more) between each pair of factors and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "cmahalanobis(mtcars, ~am + carb + index)" will print distances and plot only considering "index". Rows with NA values are omitted. 
#' 
#' @examples
#' # Example with the iris dataset
#'
#' data(iris)
#' 
#' # Calculate the Mahalanobis distance for "Species" groups in "iris" dataset
#' cmahalanobis(iris, ~Species, plot = TRUE, 
#' plot_title = "Mahalanobis Distance Between Groups", min_group_size = 3)
#'
#' # Example with the mtcars dataset
#' data(mtcars)
#' 
#' # Calculate the Mahalanobis distance for two factors in "mtcars" dataset
#' cmahalanobis(mtcars, ~am + vs, 
#' plot = TRUE, plot_title = "Mahalanobis Distance Between Groups", 
#' min_group_size = 2, pvalues_chisq = TRUE)
#' 
#' # Calculate the Mahalanobis distance for "index" in mtcars
#' cmahalanobis(mtcars, ~index, pvalues_chisq = TRUE) 
#' 
#' @export
cmahalanobis <- function(dataset, formula, plot = TRUE, 
                         plot_title = "Mahalanobis Distance Between Groups", 
                         min_group_size = 3, pvalues_chisq = FALSE) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    # Omit NA values
    dataset_imputed <- na.omit(dataset)
    
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    distan <- round(sqrt(mahalanobis(X, center = colMeans(X), cov(X))), digits = 2)
    distance_matrix <- as.matrix(distan)
    distance <- list(distance_matrix = distance_matrix)
    
    if (pvalues_chisq == TRUE) {
      p_values <- pchisq(distance_matrix^2, df = (ncol(dataset) - 1), lower.tail = FALSE)
      pval <- list(p_values = p_values)
    }
    
    if (pvalues_chisq == TRUE) {
      return(c(distance, pval))
    } else {
      return(c(distance))
    }
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  }
  else if (length(predictors) == 0) { 
    
    dataset_imputed <- na.omit(dataset)
    
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  di_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
      
      covariances <-  lapply(valid_groups, function(dt_group) {
        cov(as.matrix(dt_group[, predictors]), use = "complete.obs")
      })
      
    } else if (length(predictors) == 0) {
      
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
      
      covariances <-  lapply(valid_groups, function(dt_group) {
        cov(as.matrix(dt_group), use = "complete.obs")
      })
      
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    epsilon <- 1e-6
    for (i in 1:n) {
      mean_i <- means[[i]]
      cov_i <- covariances[[i]] + diag(epsilon, length(mean_i))
      for (j in 1:n) {
        if (i != j) {
          mean_j <- means[[j]]
          cov_j <- covariances[[j]] + diag(epsilon, length(mean_j))
          dist_ij <- round(sqrt(mahalanobis(mean_j, center = mean_i, cov = cov_i)), digits = 2)
          dist_ji <- round(sqrt(mahalanobis(mean_i, center = mean_j, cov = cov_j)), digits = 2)
          distances[i, j] <- (dist_ij + dist_ji) / 2
          distances[j,i] <- distances[i,j]
        }
      }
    }
    
    if (pvalues_chisq == TRUE) {
      d <- (distances)^2
      p_values <- round(pchisq(d, df = ((ncol(dt) - 1) * (n - 1)), lower.tail = FALSE), digits = 2)
      di_list[[grouping_var]] <- list(p_values = p_values)
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    result[[grouping_var]] <- list(distances = distances)
  }
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  if (pvalues_chisq == TRUE) {
    return(c(result, di_list))
  } else {
    return(result)
  }
}

#' @name generate_report_cmahalanobis
#' @title Generate a Microsoft Word document about the Mahalanobis distances matrix or matrices and the p-values matrix or matrices.
#' @description
#'  This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Mahalanobis distances matrix or matrices (two or more) and the p-values matrix or matrices.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Mahalanobis distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @param pvalues_chisq If TRUE, print the result of the chi-squared test on squared distances. The resulting distances with "pvalues_chisq = FALSE" are not squared; instead, with "pvalues_chisq = TRUE", the squared Mahalanobis distance matrix with corresponding p_values will be printed. Default is FALSE.
#' @return A Microsoft Word document about the Mahalanobis distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about the "Species" 
#' # factor in the iris dataset using the "permutation" method.
#' generate_report_cmahalanobis(iris, ~Species, min_group_size = 3)
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#'
#' # Generate a report about the "am" and "vs" in mtcars using "bootstrap" method.
#' generate_report_cmahalanobis(mtcars, ~am + vs,
#' pvalue.method = "bootstrap",
#' seed = 100, min_group_size = 2)
#'  
#' @export
generate_report_cmahalanobis <- function(dataset, formula, pvalue.method = "permutation",
                                         seed = NULL, min_group_size = 3, pvalues_chisq = FALSE) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if (!("index" %in% grouping_vars)) {
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    if (pvalues_chisq == FALSE) {
      cmahalanobis_results <- cmahalanobis(dataset, formula, plot = FALSE, 
                                           plot_title = "Mahalanobis Distance Between Groups", 
                                           min_group_size = min_group_size, pvalues_chisq = FALSE)
    } else if (pvalues_chisq == TRUE) {
      cmahalanobis_results <- cmahalanobis(dataset, formula, plot = FALSE, 
                                           plot_title = "Mahalanobis Distance Between Groups", 
                                           min_group_size = min_group_size, pvalues_chisq = TRUE)
    }
    
    distances <- cmahalanobis_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescmaha(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescmaha(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcmahalanobis_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cmahalanobis.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcmahalanobis.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcmahalanobis_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name pvaluescmaha
#' @title Calculate the p_values matrix or matrices (two or more) for each pair of factors inside variable or variables (two or more), using Mahalanobis distance as a base.
#' @description
#' Using the Mahalanobis distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Mahalanobis distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation" methods.
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with "airquality" dataset
#' data(airquality)
#' 
#' # Calculate p_values of "Month" variable in "airquality" dataset
#' pvaluescmaha(airquality,~Month, pvalue.method = "permutation", seed = 12,
#' min_group_size = 3)
#' 
#' # Example with "mtcars" dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" and "carb" variable in mtcars dataset
#' pvaluescmaha(mtcars,~am + carb, 
#' pvalue.method = "permutation", seed = 100, min_group_size = 2)
#' 
#' @export
pvaluescmaha <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("p_values calculation not possible. 'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("p_values calculation not possible. 'permutation' method does not have sense with 'index'")
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Mahalanobis distances using cmahalanobis
    mahalanobis_results <- cmahalanobis(dataset, as.formula(paste("~", grouping_var)),
                                        plot = FALSE, min_group_size = min_group_size)
    distances <- mahalanobis_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cmahalanobis(permuted_data,
                                         as.formula(paste("~", grouping_var)),
                                         plot = FALSE,
                                         min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cmahalanobis(bootstrap_data,
                                          as.formula(paste("~", grouping_var)),
                                          plot = FALSE,
                                          min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- mahalanobis_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
      
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name ceuclide
#' @title Calculate the Euclidean distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Euclidean distances about each pair of factors inside them. You can also select "index" to calculate the Euclidean distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Euclidean distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Euclidean distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore factors, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Euclidean distance matrix will be printed; instead, by specifying variables, the Euclidean distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "ceuclide(mtcars, ~am + carb + index)" will print distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#'
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' ceuclide(iris, ~Species, plot = TRUE, 
#' plot_title = "Euclidean Distance Between Groups", min_group_size = 2)
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' ceuclide(mtcars, ~am + carb, plot = TRUE, 
#' plot_title = "Euclidean Distance Between Groups", min_group_size = 3)
#' 
#' # Calculate ceuclide for index
#' res <- ceuclide(mtcars, ~index, 
#' min_group_size = 3)
#' 
#' @export
ceuclide <- function(dataset, formula, plot = TRUE,
                     plot_title = "Euclidean Distance Between Groups",
                     min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    # Compute the Euclidean Distance matrix based on Euclidean distance
    D <- round((dist(X)), digits = 2)
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  # List to keep track of excluded groups for each variable
  excluded_groups_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    # Identify excluded groups
    excluded_groups <- names(groups)[group_sizes < min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups (size >= ", min_group_size, ") for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    for (i in 1:n) {
      mean_i <- means[[i]]
      for (j in 1:n) {
        if (i != j) {
          mean_j <- means[[j]]
          distances[i, j] <- round(sqrt(sum2((mean_i - mean_j) ^ 2)), digits = 2)
          # The distance matrix is symmetric
          distances[j,i] <- distances[i,j]
        }
      }
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    # Save the distance matrix and excluded groups for this variable
    result[[grouping_var]] <- list(distances = distances)
  }
  
  # Print plots (if any and gridExtra is available)
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  return(result)
}

#' @name generate_report_ceuclide
#' @title Generate a Microsoft Word document about the Euclidean distances matrix or matrices and the p-values matrix or matrices.
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Euclidean distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Euclidean distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore factors, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Euclidean distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(airquality)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ceuclide(airquality, ~Month, pvalue.method = 'bootstrap',
#' min_group_size = 3)
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" and "vs" factors in mtcars dataset
#' generate_report_ceuclide(mtcars, ~am + vs, 
#' pvalue.method = 'bootstrap', seed = 100, min_group_size = 3)
#' 
#' @export
generate_report_ceuclide <- function(dataset, formula, pvalue.method = "permutation",
                                     seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    ceuclide_results <- ceuclide(dataset, formula, plot = FALSE, 
                                 plot_title = "Euclidean Distance Between Groups", 
                                 min_group_size = min_group_size)
    
    distances <- ceuclide_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluesceucl(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluesceucl(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportceuclide_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_ceuclide.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportceuclide.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportceuclide_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}


#' @name pvaluesceucl
#' @title Calculate the p_values matrix or matrices (two or more) for each pair of factors inside variable or variables (two or more), using Euclidean distance as a base.
#' @description
#' Using the Euclidean distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Euclidean distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation" methods.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesceucl(iris,~Species, pvalue.method = "permutation"
#' , min_group_size = 3)
#' 
#' # Example with mtcars dataset
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesceucl(mtcars,~am + carb, 
#' pvalue.method = "bootstrap", 
#' seed = 100, min_group_size = 2)
#' 
#' @export
pvaluesceucl <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("p_values calculation not possible. 'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("p_values calculation not possible. 'permutation' method does not have sense with 'index'")
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Euclidean distances using ceuclide
    euclidean_results <- ceuclide(dataset, as.formula(paste("~", grouping_var)),
                                  plot = FALSE, min_group_size = min_group_size)
    distances <- euclidean_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- ceuclide(permuted_data,
                                     as.formula(paste("~", grouping_var)),
                                     plot = FALSE,
                                     min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- ceuclide(bootstrap_data,
                                      as.formula(paste("~", grouping_var)),
                                      plot = FALSE,
                                      min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- euclidean_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
      
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name cmanhattan
#' @title Calculate the Manhattan distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Manhattan distances about the factors inside them. You can also select "index" to calculate the Manhattan distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Manhattan distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Manhattan distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Manhattan distances matrix will be printed; instead, by specifying variables, the Manhattan distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "cmanhattan(mtcars, ~am + carb + index)" will print the distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' cmanhattan(iris, ~Species, plot = TRUE, 
#' plot_title = "Manhattan Distance Between Groups", min_group_size = 3)
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' cmanhattan(mtcars, ~am + vs, plot = TRUE, 
#' plot_title = "Manhattan Distance Between Groups", min_group_size = 3)
#' 
#' # Calculate the Manhattan distance for 32 car models in "mtcars" dataset
#' res <- cmanhattan(mtcars, ~index, min_group_size = 3)
#' 
#' @export
cmanhattan <- function(dataset, formula, plot = TRUE, 
                       plot_title = "Manhattan Distance Between Groups", 
                       min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    # Compute the Manhattan Distance matrix based on Manhattan distance
    D <- round((dist(X, method = "manhattan")), digits = 2)
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    # Calculate Manhattan distance between each pair of groups
    for (i in 1:n) {
      mean_i <- means[[i]] # Precalculated mean of group i
      for (j in 1:n) {
        if (i != j) {
          distances[i, j] <- round(sum2(abs(mean_i - means[[j]])), digits = 2)
          distances[j,i] <- distances[i,j]
        }
      }
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    result[[grouping_var]] <- list(distances = distances)
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  return(result)
}

#' @name generate_report_cmanhattan
#' @title Generate a Microsoft Word document about the Manhattan distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Manhattan distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Manhattan distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Manhattan distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cmanhattan(iris, ~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cmanhattan(mtcars, ~am, 
#' pvalue.method = 'bootstrap', seed = 123)
#' 
#' @export
generate_report_cmanhattan <- function(dataset, formula, pvalue.method = "permutation",
                                       seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    cmanhattan_results <- cmanhattan(dataset, formula, plot = FALSE, 
                                     plot_title = "Manhattan Distance Between Groups", 
                                     min_group_size = min_group_size)
    
    distances <- cmanhattan_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescmanh(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescmanh(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcmanhattan_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cmanhattan.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcmanhattan.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcmanhattan_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name pvaluescmanh
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Manhattan distance as a base.
#' @description
#' Using the Manhattan distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Manhattan distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescmanh(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescmanh(mtcars,~am, 
#' pvalue.method = "permutation", seed = 123)
#' 
#' @export
pvaluescmanh <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Manhattan distances using cmanhattan
    manhattan_results <- cmanhattan(dataset, as.formula(paste("~", grouping_var)),
                                    plot = FALSE, min_group_size = min_group_size)
    distances <- manhattan_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cmanhattan(permuted_data,
                                       as.formula(paste("~", grouping_var)),
                                       plot = FALSE,
                                       min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cmanhattan(bootstrap_data,
                                        as.formula(paste("~", grouping_var)),
                                        plot = FALSE,
                                        min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- manhattan_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
      
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name cchebyshev
#' @title Calculate the Chebyshev distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Chebyshev distances about the factors inside them. You can also select "index" to calculate the Chebyshev distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Chebyshev distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Chebyshev distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Chebyshev distances matrix will be printed; instead, by specifying variables, the Chebyshev distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "cchebyshev(mtcars, ~am + carb + index)" will print distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples 
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' cchebyshev(iris, ~Species, plot = TRUE, 
#' plot_title = "Chebyshev Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' cchebyshev(mtcars, ~am, plot = TRUE, 
#' plot_title = "Chebyshev Distance Between Groups")
#' 
#' # Calculate the Chebyshev distance for 32 car models in "mtcars" dataset
#' res <- cchebyshev(mtcars, ~index)
#' 
#' @export
cchebyshev <- function(dataset, formula, plot = TRUE, 
                       plot_title = "Chebyshev Distance Between Groups", 
                       min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    # Compute the Chebyshev Distance matrix based on Chebyshev distance
    D <- round((dist(X, method = "maximum")), digits = 2)
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    # Calculate Chebyshev distance between each pair of groups
    for (i in 1:n) {
      mean_i <- means[[i]]  # Precalculated mean of group i
      for (j in 1:n) {
        if (i != j) {
          distances[i, j] <- round(max(abs(mean_i - means[[j]])), digits = 2)
          distances[j,i] <- distances[i,j]
        }
      }
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    result[[grouping_var]] <- list(distances = distances)
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  return(result)
}

#' @name generate_report_cchebyshev
#' @title Generate a Microsoft Word document about the Chebyshev distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Chebyshev distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Chebyshev distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Chebyshev distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cchebyshev(iris, ~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cchebyshev(mtcars, ~am, 
#' pvalue.method = "bootstrap", seed = 100)
#' 
#' @export
generate_report_cchebyshev <- function(dataset, formula, pvalue.method = "permutation",
                                       seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if (!("index" %in% grouping_vars)) {
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    cchebyshev_results <- cchebyshev(dataset, formula, plot = FALSE, 
                                     plot_title = "Chebyshev Distance Between Groups", 
                                     min_group_size = min_group_size)
    
    distances <- cchebyshev_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluesccheb(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluesccheb(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcchebyshev_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cchebyshev.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcchebyshev.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcchebyshev_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}


#' @name pvaluesccheb
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Chebyshev distance as a base.
#' @description
#' Using the Chebyshev distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Chebyshev distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".ì
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with "iris" dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccheb(iris,~Species, pvalue.method = "permutation")
#' 
#' # Example with "mtcars" dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccheb(mtcars,~am, 
#' pvalue.method = "bootstrap", seed = 100)
#' 
#' @export
pvaluesccheb <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Chebyshev distances using cchebyshev
    chebyshev_results <- cchebyshev(dataset, as.formula(paste("~", grouping_var)),
                                    plot = FALSE, min_group_size = min_group_size)
    distances <- chebyshev_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cchebyshev(permuted_data,
                                       as.formula(paste("~", grouping_var)),
                                       plot = FALSE,
                                       min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cchebyshev(bootstrap_data,
                                        as.formula(paste("~", grouping_var)),
                                        plot = FALSE,
                                        min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- chebyshev_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
      
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name chamming
#' @title Calculate the Hamming distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Hamming distances about the factors inside them. You can also select "index" to calculate the Hamming distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Hamming distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Hamming distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Hamming distances matrix will be printed; instead, by specifying variables, the Hamming distances matrix or matrices (two or more) between each pair of factors and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "chamming(mtcars, ~am + carb + index)" will print the distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' chamming(iris, ~Species, plot = TRUE, 
#' plot_title = "Hamming Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' chamming(mtcars, ~am, plot = TRUE,
#' plot_title = "Hamming Distance Between Groups")
#' 
#' # Calculate the Hamming distance for 32 car models in "mtcars" dataset
#' res <- chamming(mtcars, ~index)
#' 
#' @export
chamming <- function(dataset, formula, plot = TRUE, 
                     plot_title = "Hamming Distance Between Groups", 
                     min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, each row is its own group.  
    # Predictors are all variables except "index".
    predictor_vars <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding if necessary.
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictor_vars, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <- sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictor_vars <- colnames(predictors_data)
    
    # Compute Hamming distance for two binary vectors as the number of positions that differ.
    n_rows <- nrow(dataset_imputed)
    D <- matrix(0, nrow = n_rows, ncol = n_rows)
    for (i in 1:(n_rows - 1)) {
      for (j in (i + 1):n_rows) {
        hamming_dist <- sum2(predictors_data[i, ] != predictors_data[j, ])
        D[i, j] <- hamming_dist
        D[j, i] <- hamming_dist
      }
    }
    rownames(D) <- colnames(D) <- as.character(dataset_imputed$index)
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
    
  }
  
  # Otherwise, Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n_groups <- length(valid_groups)
    
    # Compute pairwise Hamming distances between the group binary vectors.
    distances <- matrix(0, nrow = n_groups, ncol = n_groups)
    rownames(distances) <- colnames(distances) <- group_names
    
    for (i in 1:(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        a <- valid_groups[[j]]
        b <- valid_groups[[i]]
        if ((nrow(a)) < (nrow(b))) {
          b <- b[1:nrow(a), ]
          hamming_dist <- sum2(b != a)
        } else if ((nrow(a) > (nrow(b)))) {
          a <- a[1:nrow(b), ]
          hamming_dist <- sum2(b != a)
        } else if ((nrow(b)) == (nrow(a))) {
          hamming_dist <- sum2(b != a)
        }
        distances[i, j] <- hamming_dist 
        distances[j, i] <- distances[i, j]
      }
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    result[[grouping_var]] <- list(distances = distances)
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  return(result)
}


#' @name pvalueschamm
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Hamming distance as a base.
#' @description
#' Using the Hamming distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Hamming distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with "iris" dataset
#' 
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in "iris" dataset
#' pvalueschamm(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#'  
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvalueschamm(mtcars,~am, 
#' pvalue.method = "permutation", seed = 100)
#' 
#' @export
pvalueschamm <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Hamming distances using chamming
    hamming_results <- chamming(dataset, as.formula(paste("~", grouping_var)),
                                plot = FALSE, min_group_size = min_group_size)
    distances <- hamming_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- chamming(permuted_data,
                                     as.formula(paste("~", grouping_var)),
                                     plot = FALSE,
                                     min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- chamming(bootstrap_data,
                                      as.formula(paste("~", grouping_var)),
                                      plot = FALSE,
                                      min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- hamming_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
      
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}


#' @name generate_report_chamming
#' @title Generate a Microsoft Word document about the Hamming distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Hamming distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Hamming distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Hamming distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_chamming(iris, ~Species)
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_chamming(mtcars, ~am, 
#' pvalue.method = "bootstrap", seed = 124)
#' 
#' @export
generate_report_chamming <- function(dataset, formula, pvalue.method = "permutation",
                                     seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    chamming_results <- chamming(dataset, formula, plot = FALSE, 
                                 plot_title = "Hamming Distance Between Groups", 
                                 min_group_size = min_group_size)
    
    distances <- chamming_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvalueschamm(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvalueschamm(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportchamming_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_chamming.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportchamming.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportchamming_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name ccanberra
#' @title Calculate the Canberra distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Canberra distances about the factors inside them. You can also select "index" to calculate the Canberra distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Canberra distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Canberra distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Canberra distances matrix will be printed; instead, by specifying variables, the Canberra distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "ccanberra(mtcars, ~am + carb + index)" will print distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' 
#' data(iris)
#'
#' ccanberra(iris, ~Species, plot = TRUE, 
#' plot_title = "Canberra Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' 
#' data(mtcars)
#' 
#' ccanberra(mtcars, ~am, plot = TRUE, 
#' plot_title = "Canberra Distance Between Groups")
#' 
#' # Calculate the Canberra distance for 32 car models in "mtcars" dataset
#' res <- ccanberra(mtcars, ~index)
#' 
#' @export
ccanberra <- function(dataset, formula, plot = TRUE, 
                      plot_title = "Canberra Distance Between Groups", 
                      min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    # Compute the Canberra Distance matrix based on standard Canberra distance
    D <- as.matrix(round(dist(X, method = "canberra"), digits = 2))
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    # Calculate Canberra distance between each pair of groups.
    for (i in 1:n) {
      mean_i <- means[[i]]
      for (j in 1:n) {
        if (i != j) {
          mean_j <- means[[j]]
          mean_diff <- abs(mean_i - mean_j)
          mean_sum <- abs(mean_i) + abs(mean_j)
          distances[i,j] <- round(sum2(mean_diff / mean_sum), digits = 2)
          distances[j, i] <- distances[i, j]  # Ensure symmetry
        }
      }
    }
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    result[[grouping_var]] <- list(distances = distances)
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  return(result)
}

#' @name pvaluesccanb
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Canberra distance as a base.
#' @description
#' Using the Canberra distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Canberra distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccanb(iris,~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccanb(mtcars,~am + vs, 
#' pvalue.method = "permutation", seed = 100)
#' 
#' @export
pvaluesccanb <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Canberra distances using ccanberra
    canberra_results <- ccanberra(dataset, as.formula(paste("~", grouping_var)),
                                  plot = FALSE, min_group_size = min_group_size)
    distances <- canberra_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- ccanberra(permuted_data,
                                      as.formula(paste("~", grouping_var)),
                                      plot = FALSE,
                                      min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- ccanberra(bootstrap_data,
                                       as.formula(paste("~", grouping_var)),
                                       plot = FALSE,
                                       min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- canberra_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}


#' @name generate_report_ccanberra
#' @title Generate a Microsoft Word document about the Canberra distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Canberra distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Canberra distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Canberra distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ccanberra(iris, ~Species, 
#' pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_ccanberra(mtcars, ~am, 
#' pvalue.method = "bootstrap", seed = 123)
#' 
#' # Generate a report for 32 car models in "mtcars" dataset,
#' # using "bootstrap" method
#' generate_report_ccanberra(mtcars, ~am, pvalue.method = "bootstrap")
#' 
#' @export
generate_report_ccanberra <- function(dataset, formula, pvalue.method = "permutation",
                                      seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    ccanberra_results <- ccanberra(dataset, formula, plot = FALSE, 
                                   plot_title = "Canberra Distance Between Groups", 
                                   min_group_size = min_group_size)
    
    distances <- ccanberra_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluesccanb(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluesccanb(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportccanberra_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_ccanberra.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportccanberra.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportccanberra_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name cminkowski
#' @title Calculate the Minkowski distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Minkowski distances about the factors inside them. You can also select "index" to calculate the Minkowski distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Minkowski distances matrix or matrices (two or more).
#' @param p Order of the Minkowski distance.
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Minkowski distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Minkowski distances matrix will be printed; instead, by specifying variables, the Minkowski distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' When p < 1, the Minkowski distance is a "dissimilarity" measure. When p >= 1, 
#' the triangle inequality property is satisfied and we say "Minkowski distance". If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "cminkowski(mtcars, ~am + carb + index)" will print distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' cminkowski(iris, ~Species, p = 3, plot = TRUE, 
#' plot_title = "Minkowski Distance Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' cminkowski(mtcars, ~am, p = 3, plot = TRUE, 
#' plot_title = "Minkowski Distance Between Groups")
#' 
#' # Calculate the Minkowski distance for 32 car models in "mtcars" dataset
#' res <- cminkowski(mtcars, ~index, p = 2, plot = TRUE)
#' 
#' @export
cminkowski <- function(dataset, formula, p = 3, plot = TRUE, 
                       plot_title = "Minkowski Distance Between Groups", 
                       min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    # Compute the Minkowski dissimilarity/distance matrix based on standard Minkowski distance
    D <- as.matrix(round(dist(X, method = "minkowski", p = p), digits = 2))
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    # Calculate Minkowski dissimilarity/distance between each pair of groups.
    for (i in 1:n) {
      mean_i <- means[[i]]
      for (j in 1:n) {
        if (i != j) {
          mean_j <- means[[j]]
          # The Minkowski dissimilarity/distance
          distances[i, j] <- round( sum2(abs(mean_i - mean_j)^(p)) ^(1/p), digits = 2)
          distances[j, i] <- distances[i, j]  # Ensure symmetry
        }
      }
    }
    
    if (p >= 1) {
      
      if (plot) {
        plot_title_i <- paste(plot_title, "-", grouping_var)
        
        dist_df <- as.matrix((distances))
        dist_df1 <- reshape2::melt(dist_df)
        
        f <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
          ggplot2::geom_tile(color = "white") +
          ggplot2::scale_colour_brewer(palette = "Greens") +
          ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
          ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
          ggplot2::theme_minimal()
        
        plot_list[[grouping_var]] <- f
      }
      
      result[[grouping_var]] <- list(distances = distances)
      if (plot && length(plot_list) > 0) {
        gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
      }
    }
    else if (p < 1) {
      if (plot) {
        plot_title <- "Minkowski Dissimilarity Between Groups"
        
        dist_df <- as.matrix((distances))
        dist_df1 <- reshape2::melt(dist_df)
        
        f <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
          ggplot2::geom_tile(color = "white") +
          ggplot2::scale_colour_brewer(palette = "Greens") +
          ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
          ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
          ggplot2::theme_minimal()
        
        
        plot_list[[grouping_var]] <- f
      }
      
      result[[grouping_var]] <- list(distances = distances)
      if (plot && length(plot_list) > 0) {
        gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
      }
    }
  }
  return(result)
}

#' @name pvaluescmink
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Minkowski dissimilarity/distance as a base.
#' @description
#' Using the Minkowski dissimilarity/distance for the dissimilarities/distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Minkowski dissimilarities/distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param p Order of the Minkowski dissimilarities/distances. The default value is 3.
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' 
#' # data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' 
#' pvaluescmink(iris,~Species, p = 3, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' 
#' pvaluescmink(mtcars,~am, p = 3, 
#' pvalue.method = "permutation", seed = 100)
#' 
#' @note
#' When p < 1, the Minkowski distance is a "dissimilarity" measure. When p >= 1, 
#' the triangle inequality property is satisfied and we say "Minkowski distance".
#' 
#' @export
pvaluescmink <- function(dataset, formula, pvalue.method = "permutation", p = 3, plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Minkowski distances/dissimilarities using cminkowski
    minkowski_results <- cminkowski(dataset, as.formula(paste("~", grouping_var)), p = p,
                                    plot = FALSE, min_group_size = min_group_size)
    distances <- minkowski_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cminkowski(permuted_data,
                                       as.formula(paste("~", grouping_var)), p = p,
                                       plot = FALSE,
                                       min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cminkowski(bootstrap_data,
                                        as.formula(paste("~", grouping_var)), p = p,
                                        plot = FALSE,
                                        min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- minkowski_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_cminkowski
#' @title Generate a Microsoft Word document about the Minkowski dissimilarities/distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Minkowski dissimilarities/distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Minkowski dissimilarities/distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param p Order of the Minkowski dissimilarities/distances. The default value is 3.
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Minkowski dissimilarities/distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cminkowski(iris, ~Species, p = 3, 
#' pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cminkowski(mtcars, ~am, 
#' p = 3, pvalue.method = 'permutation', seed = 234)
#'  
#' @details
#' When p < 1, the Minkowski distance is a "dissimilarity" measure. When p >= 1, 
#' the triangle inequality property is satisfied and we say "Minkowski distance".
#' 
#' @export
generate_report_cminkowski <- function(dataset, formula, p = 3, pvalue.method = "permutation",
                                       seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    cminkowski_results <- cminkowski(dataset, formula, p = p, plot = FALSE, 
                                     plot_title = "Minkowski Distance Between Groups", 
                                     min_group_size = min_group_size)
    
    distances <- cminkowski_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescmink(dataset, formula, p = p, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescmink(dataset, formula, p = p, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcminkowski_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cminkowski.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcminkowski.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcminkowski_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name ccosine
#' @title Calculate the Cosine dissimilarities for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Cosine dissimilarities about the factors inside them. You can also select "index" to calculate the Cosine dissimilarities between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Cosine dissimilarities matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Cosine dissimilarities matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Cosine dissimilarities matrix will be printed; instead, by specifying variables, the Cosine dissimilarities matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only dissimilarities between rows are calculated. Therefore, this snippet: "ccosine(mtcars, ~am + carb + index)" will print dissimilarities only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with iris dataset
#' 
#' data(iris)
#' 
#' ccosine(iris, ~Species, plot = TRUE, 
#' plot_title = "Cosine Dissimilarity Between Groups")
#' 
#' # Example with mtcars dataset
#' 
#' data(mtcars)
#' 
#' ccosine(mtcars, ~am, plot = TRUE, 
#' plot_title = "Cosine Dissimilarity Between Groups")
#' 
#' # Calculate the Cosine dissimilarity for 32 car models in "mtcars" dataset
#' res <- ccosine(mtcars, ~index)
#' 
#' @export
ccosine <- function(dataset, formula, plot = TRUE, 
                    plot_title = "Cosine Dissimilarity Between Groups", 
                    min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # If the user specified "~index", use the individual mode.
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # In this case, grouping is performed for each row (i.e., each observation is its own group).
    # Predictors are all variables except for "index".
    predictors <- setdiff(names(dataset), "index")
    
    ## Perform imputation and one-hot encoding (if necessary)
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Convert predictors to a matrix to speed up calculations and standardize
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    X <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    n_rows <- nrow(X)
    D <- matrix(0, nrow = n_rows, ncol = n_rows)
    for (i in (1:n_rows)) {
      for (j in (i:n_rows)) {
        similarity <- crossprod(X[i,], X[j,]) / sqrt((crossprod(X[i,])) * (crossprod(X[j,])))
        dissimilarity <- 1 - similarity
        D[i,j] <- round(dissimilarity, digits = 2)
        D[j,i] <- D[i,j]
      }
    }
    rownames(D) <- colnames(D) <- as.character(dataset_imputed$index)
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    
    # Calculate the cosine dissimilarity between each pair of groups' mean vectors
    for (i in 1:n) {
      for (j in i:n) {
        if (i != j) {
          similarity <- crossprod(means[[i]], means[[j]]) / sqrt((crossprod(means[[i]])) * (crossprod(means[[j]])))
          distances[i, j] <- round(1 - similarity, digits = 2)
          distances[j, i] <- distances[i, j]
        }
      }
    }
    
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p
    }
    
    result[[grouping_var]] <- list(distances = distances)
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  
  return(result)
}

#' @name pvaluesccosi
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Cosine dissimilarity as a base.
#' @description
#' Using the Cosine dissimilarity for the dissimilarities calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Cosine dissimilarities matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size  Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluesccosi(iris,~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluesccosi(mtcars,~am, 
#' pvalue.method = "permutation", seed = 123)
#' 
#' @export
pvaluesccosi <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Cosine dissimilarity using ccosine
    cosine_results <- ccosine(dataset, as.formula(paste("~", grouping_var)),
                              plot = FALSE, min_group_size = min_group_size)
    distances <- cosine_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- ccosine(permuted_data,
                                    as.formula(paste("~", grouping_var)),
                                    plot = FALSE,
                                    min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- ccosine(bootstrap_data,
                                     as.formula(paste("~", grouping_var)),
                                     plot = FALSE,
                                     min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- cosine_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_ccosine
#' @title Generate a Microsoft Word document about the Cosine dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Cosine dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Cosine dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Cosine dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_ccosine(iris, ~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_ccosine(mtcars, ~am,
#' pvalue.method = "bootstrap", seed = 123)
#' 
#' @export
generate_report_ccosine <- function(dataset, formula, pvalue.method = "permutation",
                                    seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    ccosine_results <- ccosine(dataset, formula, plot = FALSE, 
                               plot_title = "Cosine Dissimilarity Between Groups", 
                               min_group_size = min_group_size)
    
    distances <- ccosine_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluesccosi(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluesccosi(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportccosine_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_ccosine.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportccosine.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportccosine_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name cbhattacharyya
#' @title Calculate the Bhattacharyya dissimilarities for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Bhattacharyya dissimilarities about the factors inside them. You can also select "index" to calculate the Bhattacharyya dissimilarities between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Bhattacharyya dissimilarities matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Bhattacharyya dissimilarities matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Bhattacharyya dissimilarities matrix will be printed; instead, by specifying variables, the Bhattacharyya dissimilarities matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only dissimilarities between rows are calculated. Therefore, this snippet: "cbhattacharyya(mtcars, ~am + carb + index)" will print dissimilarities only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' 
#' data(iris)
#' 
#' cbhattacharyya(iris, ~Species, plot = TRUE, 
#' plot_title = "Bhattacharyya Dissimilarity Between Groups")
#'
#' # Example with the mtcars dataset
#' 
#' data(mtcars)
#' 
#' cbhattacharyya(mtcars, ~am, plot = TRUE, 
#' plot_title = "Bhattacharyya Dissimilarity Between Groups")
#' 
#' # Calculate Bhattacharyya distance for index
#' res <- cbhattacharyya(mtcars, ~index)
#' 
#' @export
cbhattacharyya <- function(dataset, formula, plot = TRUE,
                           plot_title = "Bhattacharyya Dissimilarity Between Groups",
                           min_group_size = 3) {
  
  if (!is.data.frame(dataset)) stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) stop("At least one grouping variable must be specified in the formula.")
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) stop("Some grouping variables are not present in the dataset.")
  
  
  predictor_vars <- setdiff(names(dataset), grouping_vars)
  
  dataset_imputed <- na.omit(dataset)
  
  predictors_data <- dataset_imputed[, predictor_vars, drop = FALSE]
  predictors_data <- model.matrix(~ . - 1, data = predictors_data)
  predictors_data <- as.data.frame(predictors_data)
  
  variances <- sapply(predictors_data, stats::var, na.rm = TRUE)
  predictors_data <- predictors_data[, !is.na(variances) & variances > 1e-9, drop = FALSE]
  predictor_vars <- colnames(predictors_data)
  
  # Index mode
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    data_matrix <- as.matrix(predictors_data)
    n_rows <- nrow(data_matrix)
    D <- matrix(0, nrow = n_rows, ncol = n_rows)
    rownames(D) <- colnames(D) <- dataset$index
    
    # Normalize each row to sum to 1
    row_sums <- rowSums2(data_matrix, na.rm = TRUE)
    # Avoid division by zero if a row is completely zero
    row_sums[row_sums == 0] <- 1
    distributions <- data_matrix / row_sums
    
    for (i in 1:n_rows) {
      for (j in (i:n_rows)) {
        if (i == j) next
        
        p <- distributions[i, ]
        q <- distributions[j, ]
        
        bhattacharyya_coeff <- sum2((sqrt(p * q)))
        bhattacharyya_dist <- -log(bhattacharyya_coeff)
        
        D[i, j] <- round(bhattacharyya_dist, digits = 2)
        D[j, i] <- D[i, j]
      }
    }
    
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  # Group mode
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) { 
    dataset_prepared <- cbind(dataset_imputed[, grouping_vars, drop = FALSE], predictors_data)
    dt <- as.data.frame(dataset_prepared)
  } else if (length(predictors) == 0) {
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (g in grouping_vars) {
    groups <- split(dt, dt[[g]])
    group_sizes <- sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", g, "- skipping."))
      next
    }
    group_names <- names(valid_groups)
    
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        # Normalize each row to sum to 1
        row_sums <- rowSums2(as.matrix(dt_group[, predictors]))
        # Avoid division by zero if a row is completely zero
        row_sums[row_sums == 0] <- 1
        distributions <- (dt_group[, predictors]) / row_sums
        colMeans2(as.matrix(distributions))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        # Normalize each row to sum to 1
        row_sums <- rowSums2(as.matrix(dt_group))
        # Avoid division by zero if a row is completely zero
        row_sums[row_sums == 0] <- 1
        distributions <- dt_group / row_sums
        colMeans2(as.matrix(distributions))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        mean_i <- means[[i]]
        mean_j <- means[[j]]
        
        bhattacharyya_coeff <- sum2((sqrt(mean_i * mean_j)))
        bhattacharyya_dist <- -log(bhattacharyya_coeff)
        
        distances[i, j] <- round(bhattacharyya_dist, digits = 2)
        distances[j, i] <- distances[i, j]
      }
    }
    
    result[[g]] <- list(distances = distances)
    if (plot) {
      plot_title_g <- paste(plot_title, "-", g)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p_obj <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[g]] <- p_obj
    }
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  
  return(result)
}

#' @name pvaluescbatt
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Bhattacharyya dissimilarities as a base.
#' @description
#' Using the Bhattacharyya dissimilarity for the dissimilarities calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate Bhattacharyya dissimilarities matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescbatt(iris,~Species, pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#'
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescbatt(mtcars,~am, 
#' pvalue.method = "bootstrap", seed = 123)
#' 
#' @export
pvaluescbatt <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Bhattacharyya distances using cbhattacharyya
    bhattacharyya_results <- cbhattacharyya(dataset, as.formula(paste("~", grouping_var)),
                                            plot = FALSE, min_group_size = min_group_size)
    distances <- bhattacharyya_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cbhattacharyya(permuted_data,
                                           as.formula(paste("~", grouping_var)),
                                           plot = FALSE,
                                           min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cbhattacharyya(bootstrap_data,
                                            as.formula(paste("~", grouping_var)),
                                            plot = FALSE,
                                            min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- bhattacharyya_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_cbhattacharyya
#' @title Generate a Microsoft Word document about the Bhattacharyya dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Bhattacharyya dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Bhattacharyya dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size  Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Bhattacharyya dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cbhattacharyya(iris, ~Species, 
#' pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cbhattacharyya(mtcars, ~am, 
#' pvalue.method = "bootstrap", seed = 123)
#' 
#' @export
generate_report_cbhattacharyya <- function(dataset, formula, pvalue.method = "permutation",
                                           seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    bhattacharyya_results <- cbhattacharyya(dataset, formula, plot = FALSE, 
                                            plot_title = "Bhattacharyya Dissimilarity Between Groups", 
                                            min_group_size = min_group_size)
    
    distances <- bhattacharyya_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescbatt(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescbatt(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcbhattacharyya_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cbhattacharyya.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcbhattacharyya.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportccbhattacharyya_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name cjaccard
#' @title Calculate the Jaccard distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Jaccard distances about the factors inside them. You can also select "index" to calculate the Jaccard distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Jaccard distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Jaccard distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Jaccard distances matrix will be printed; instead, by specifying variables, the Jaccard distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "cjaccard(mtcars, ~am + carb + index)" will print distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' 
#' data(iris)
#'
#' cjaccard(iris, ~Species, plot = TRUE,
#' plot_title = "Jaccard Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' 
#' data(mtcars)
#' 
#' cjaccard(mtcars, ~am, 
#' plot = TRUE, plot_title = "Jaccard Distance Between Groups")
#' 
#' res <- cjaccard(mtcars, ~index,
#' plot = TRUE)
#' 
#' @export
cjaccard <- function(dataset, formula, plot = TRUE,
                     plot_title = "Jaccard Distance Between Groups",
                     min_group_size = 3) {
  
  if (!is.data.frame(dataset)) stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars)) stop("Some grouping variables are not present in the dataset.")
  
  # Index mode
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    predictor_vars <- setdiff(names(dataset), "index")
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictor_vars, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <- sapply(predictors_data, stats::var, na.rm = TRUE)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    # Compute the Jaccard Distance matrix based on standard Jaccard distance
    D <-  as.matrix(round(dist(predictors_data, method = "binary"), digits = 2))
    
    return(list(index = list(distances = D)))
  }
  
  # Group mode
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (g in grouping_vars) {
    groups <- split(dt, dt[[g]])
    group_sizes <- sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", g, "- skipping."))
      next
    }
    group_names <- names(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    n <- length(valid_groups)
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    
    for (i in 1:n) {
      for (j in i:n) { 
        if (i == j) next 
        
        mean_i <- means[[i]]
        mean_j <- means[[j]]
        
        intersection_val <- (length(intersect(mean_i, mean_j)))
        union_val <- (length(union(mean_i, mean_j)))
        
        jaccard_similarity <- intersection_val / union_val
        jaccard_distance <- 1 - jaccard_similarity
        
        distances[i, j] <- round(jaccard_distance, digits = 2)
        distances[j, i] <- distances[i, j]
      }
    }
    
    result[[g]] <- list(distances = distances)
    if (plot) {
      plot_title_g <- paste(plot_title, "-", g)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p_obj <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[g]] <- p_obj
    }
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  
  return(result)
}

#' @name pvaluescjacc
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Jaccard distance as a base.
#' @description
#' Using the Jaccard distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Jaccard distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with the iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescjacc(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with the mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescjacc(mtcars,~am, 
#' pvalue.method = "permutation",
#'  seed = 122)
#' 
#' @export
pvaluescjacc <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Jaccard distances using cjaccard
    jaccard_results <- cjaccard(dataset, as.formula(paste("~", grouping_var)),
                                plot = FALSE, min_group_size = min_group_size)
    distances <- jaccard_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cjaccard(permuted_data,
                                     as.formula(paste("~", grouping_var)),
                                     plot = FALSE,
                                     min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cjaccard(bootstrap_data,
                                      as.formula(paste("~", grouping_var)),
                                      plot = FALSE,
                                      min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- jaccard_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_cjaccard
#' @title Generate a Microsoft Word document about the Jaccard distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor and returns a Microsoft Word document about the Jaccard distances matrix or matrices and the p-values matrix or matrices.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Jaccard distances matrix or matrices and the p_values matrix or matrices.
#' @param pvalue.method A p_value method used to calculate the matrix, the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Jaccard distance matrix or matrices and the p_values matrix or matrices.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cjaccard(iris, ~Species,
#' pvalue.method = "permutation")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cjaccard(mtcars, ~am,
#' pvalue.method = "bootstrap", seed = 223)
#' 
#' @export
generate_report_cjaccard <- function(dataset, formula, pvalue.method = "permutation",
                                     seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    jaccard_results <- cjaccard(dataset, formula, plot = FALSE, 
                                plot_title = "Jaccard Distance Between Groups", 
                                min_group_size = min_group_size)
    
    distances <- jaccard_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescjacc(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescjacc(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcjaccard_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cjaccard.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcjaccard.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcjaccard_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name chellinger
#' @title Calculate the Hellinger distances for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Hellinger distances about the factors inside them. You can also select "index" to calculate the Hellinger distances between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Hellinger distances matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Hellinger distances matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Hellinger distances matrix will be printed; instead, by specifying variables, the Hellinger distances matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only distances between rows are calculated. Therefore, this snippet: "chellinger(mtcars, ~am + carb + index)" will print the distances only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' 
#' data(iris)
#'
#' chellinger(iris, ~Species, plot = TRUE,
#' plot_title = "Hellinger Distance Between Groups")
#'
#' # Example with the mtcars dataset
#' 
#' data(mtcars)
#' 
#' chellinger(mtcars, ~am, plot = TRUE,
#' plot_title = "Hellinger Distance Between Groups")
#' 
#' res <- chellinger(mtcars, ~index)
#' 
#' @export
chellinger <- function(dataset, formula, plot = TRUE,
                       plot_title = "Hellinger Distance Between Groups",
                       min_group_size = 3) {
  
  if (!is.data.frame(dataset)) stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars)) stop("Some grouping variables are not present in the dataset.")
  
  predictor_vars <- setdiff(names(dataset), grouping_vars)
  
  dataset_imputed <- na.omit(dataset)
  predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictor_vars, drop = FALSE])
  predictors_data <- as.data.frame(predictors_data)
  
  variances <- sapply(predictors_data, stats::var, na.rm = TRUE)
  predictors_data <- predictors_data[, !is.na(variances) & variances > 1e-9, drop = FALSE]
  predictor_vars <- colnames(predictors_data)
  
  # Index mode
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    # Prepare the predictors matrix
    data_matrix <- as.matrix(predictors_data)
    n_rows <- nrow(data_matrix)
    D <- matrix(0, nrow = n_rows, ncol = n_rows)
    
    # Normalize each row to sum to 1
    row_sums <- rowSums2(data_matrix)
    # Avoid division by zero if a row is completely zero
    row_sums[row_sums == 0] <- 1
    distributions <- data_matrix / row_sums
    
    rownames(D) <- colnames(D) <- seq(from = 1, by = 1, length.out = nrow(dataset_imputed))
    
    for (i in 1:(n_rows)) {
      for (j in (i:n_rows)) {
        
        i1 <- distributions[i,]
        i2 <- distributions[j,]
        
        h <- ((1/sqrt(2)) * sqrt(sum((sqrt(i1) - sqrt(i2))^2)))
        D[i, j] <- round(h, digits = 2)
        D[j, i] <- D[i,j]
      }
    }
    
    D <- as.matrix(D)
    
    return(list(index = list(distances = D)))
  }
  
  # Group mode
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    dataset_prepared <- cbind(dataset_imputed[, grouping_vars, drop = FALSE], predictors_data)
    dt <- as.data.frame(dataset_prepared)
  } else if (length(predictors) == 0) {
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  
  for (g in grouping_vars) {
    groups <- split(dt, dt[[g]])
    group_sizes <- sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", g, "- skipping."))
      next
    }
    group_names <- names(valid_groups)
    
    n <- length(valid_groups)
    distances <- matrix(0, nrow = n, ncol = n)
    rownames(distances) <- colnames(distances) <- group_names
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        # Normalize each row to sum to 1
        row_sums <- rowSums2(as.matrix(dt_group[, predictors]))
        # Avoid division by zero if a row is completely zero
        row_sums[row_sums == 0] <- 1
        distributions <- (dt_group[, predictors]) / row_sums
        colMeans2(as.matrix(distributions))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        # Normalize each row to sum to 1
        row_sums <- rowSums2(as.matrix(dt_group))
        # Avoid division by zero if a row is completely zero
        row_sums[row_sums == 0] <- 1
        distributions <- dt_group / row_sums
        colMeans2(as.matrix(distributions))
      })
    }
    
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        p_group <- means[[i]]
        q_group <- means[[j]]
        
        # Hellinger calculation
        sum_sqrt_group <- sqrt(sum((sqrt(p_group) - sqrt(q_group))^2))
        hellinger_dist_group <- sum_sqrt_group * (1/(sqrt(2)))
        
        distances[i, j] <- round(hellinger_dist_group, digits = 2)
        distances[j, i] <- distances[i, j]
      }
    }
    
    result[[g]] <- list(distances = distances)
    if (plot) {
      plot_title_g <- paste(plot_title, "-", g)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p_obj <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[g]] <- p_obj
    }
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  
  return(result)
}

#' @name pvalueschell
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Hellinger distances as a base.
#' @description
#' Using the Hellinger distance for the distances calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Hellinger distances matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for 'bootstrap' and 'permutation'.
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvalueschell(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvalueschell(mtcars,~am,
#' pvalue.method = "permutation", seed = 122)
#' 
#' @export
pvalueschell <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Hellinger distances using chellinger
    hellinger_results <- chellinger(dataset, as.formula(paste("~", grouping_var)),
                                    plot = FALSE, min_group_size = min_group_size)
    distances <- hellinger_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- chellinger(permuted_data,
                                       as.formula(paste("~", grouping_var)),
                                       plot = FALSE,
                                       min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- chellinger(bootstrap_data,
                                        as.formula(paste("~", grouping_var)),
                                        plot = FALSE,
                                        min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- hellinger_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_chellinger
#' @title Generate a Microsoft Word document about the Hellinger distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Hellinger distances matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Hellinger distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for 'bootstrap' and 'permutation'.
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Hellinger distances matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_chellinger(iris, ~Species,
#' pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_chellinger(mtcars, ~am, 
#' pvalue.method = "bootstrap", seed = 100)
#' 
#' @export
generate_report_chellinger <- function(dataset, formula, pvalue.method = "permutation",
                                       seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    hellinger_results <- chellinger(dataset, formula, plot = FALSE, 
                                    plot_title = "Hellinger Distance Between Groups", 
                                    min_group_size = min_group_size)
    
    distances <- hellinger_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvalueschell(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvalueschell(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportchellinger_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_chellinger.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportchellinger.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportchellinger_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name cbraycurtis
#' @title Calculate the Bray-Curtis dissimilarities for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Bray-Curtis dissimilarities about the factors inside them. You can also select "index" to calculate the Bray-Curtis dissimilarities between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Bray-Curtis dissimilarities matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Bray-Curtis dissimilarities matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Bray-Curtis dissimilarities matrix will be printed; instead, by specifying variables, the Bray-Curtis dissimilarities matrix or matrices (two or more) between each pair of factors and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only dissimilarities between rows are calculated. Therefore, this snippet: "cbraycurtis(mtcars, ~am + carb + index)" will print dissimilarities only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' data(iris)
#'
#' cbraycurtis(iris, ~Species, plot = TRUE,
#' plot_title = "Bray-Curtis Dissimilarity Between Groups")
#'
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Example with the mtcars dataset
#' cbraycurtis(mtcars, ~am, 
#' plot = TRUE, plot_title = "Bray-Curtis Dissimilarity Between Groups")
#' 
#' # Calculate the Bray-Curtis dissimilarity for 32 car models in "mtcars" dataset
#' res <- cbraycurtis(mtcars, ~index)
#' 
#' @export
cbraycurtis <- function(dataset, formula, plot = TRUE, 
                        plot_title = "Bray-Curtis Dissimilarity Between Groups", 
                        min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # Index mode
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    predictors <- setdiff(names(dataset), "index")
    
    dataset_imputed <- na.omit(dataset)
    
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    
    n <- nrow(predictors_data)
    
    D <- matrix(0, nrow = n, ncol = n)
    
    for (i in (1:n)) {
      for (j in (i:n)) {
        sim <- (sum(abs(predictors_data[i,] - predictors_data[j,]))) / sum((predictors_data[i,] + predictors_data[j,]))
        D[i, j] <- round(sim, digits = 2)
        D[j, i] <- D[i, j]
      }
    }
    
    D <- as.matrix(D)
    
    rownames(D) <- colnames(D) <- as.character(dataset_imputed$index)
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    
    n <- length(valid_groups)
    distances <- matrix(0, nrow = n, ncol = n)
    
    rownames(distances) <- colnames(distances) <- group_names
    
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        vec_i <- means[[i]]
        vec_j <- means[[j]]
        
        similarity <- (sum(abs(vec_i - vec_j))) / sum((vec_i + vec_j))
        
        distances[i, j] <- round(similarity, digits = 2)
        distances[j, i] <- distances[i, j]
      }
    }
    
    result[[grouping_var]] <- list(distances = distances)
    if (plot) {
      plot_title_g <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p_obj <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p_obj
    }
  }
  
  if (plot && length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ceiling(sqrt(length(plot_list))))
  }
  
  return(result)
}

#' @name pvaluescbrcu
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Bray-Curtis dissimilarity as a base.
#' @description
#' Using the Bray-Curtis dissimilarity for the dissimilarities calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Bray-Curtis dissimilarities matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for 'bootstrap' or 'permutation'.
#' @param min_group_size Minimum group size to maintain. The default value is 3,therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescbrcu(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescbrcu(mtcars,~am,
#' pvalue.method = "permutation", seed = 111)
#' 
#' @export
pvaluescbrcu <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Bray-Curtis dissimilarities using cbraycurtis
    braycurtis_results <- cbraycurtis(dataset, as.formula(paste("~", grouping_var)),
                                      plot = FALSE, min_group_size = min_group_size)
    distances <- braycurtis_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- cbraycurtis(permuted_data,
                                        as.formula(paste("~", grouping_var)),
                                        plot = FALSE,
                                        min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- cbraycurtis(bootstrap_data,
                                         as.formula(paste("~", grouping_var)),
                                         plot = FALSE,
                                         min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean2(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- braycurtis_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_cbraycurtis
#' @title Generate a Microsoft Word document about the Bray-Curtis dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Bray-Curtis dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Bray-Curtis dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for 'bootstrap' or 'permutation'.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Bray-Curtis dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_cbraycurtis(iris, ~Species, 
#' pvalue.method = "permutations")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_cbraycurtis(mtcars, ~am, 
#' pvalue.method = 'bootstrap', seed = 124)
#' 
#' @export
generate_report_cbraycurtis <- function(dataset, formula, pvalue.method = "permutation",
                                        seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  
  if (!("index" %in% grouping_vars)) {
    
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    braycurtis_results <- cbraycurtis(dataset, formula, plot = FALSE, 
                                      plot_title = "Bray-Curtis Dissimilarity Between Groups", 
                                      min_group_size = min_group_size)
    
    distances <- braycurtis_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescbrcu(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescbrcu(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcbraycurtis_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_cbraycurtis.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcbraycurtis.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcbraycurtis_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}

#' @name csorensendice
#' @title Calculate the Sorensen-Dice dissimilarities for each pair of factors or for the index.
#' @description
#' This function takes a dataframe and a variable or variables (two or more) in input, and returns a matrix or matrices (two or more) with the Sorensen-Dice dissimilarities about the factors inside them. You can also select "index" to calculate the Sorensen-Dice dissimilarities between each row. 
#' @param dataset A dataframe.
#' @param formula The index of the dataframe, otherwise a variable or variables (two or more) with factors which you want to calculate the Sorensen-Dice dissimilarities matrix or matrices (two or more). 
#' @param plot Logical, if TRUE, a plot or plots (two or more) of the Sorensen-Dice dissimilarities matrix or matrices about factors (two or more) are displayed.
#' @param plot_title If plot is TRUE, the title to be used for plot or plots about factors. The default value is TRUE.
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded. For "index", this value is always 1.
#' @return According to the option chosen in formula, with "index" the Sorensen-Dice dissimilarities matrix will be printed; instead, by specifying variables, the Sorensen-Dice dissimilarities matrix or matrices (two or more) between each pair of groups and, optionally, the plot or plots (two or more) will be printed.
#' 
#' @note
#' If "index" is selected with variables, only dissimilarities between rows are calculated. Therefore, this snippet: "csorensendice(mtcars, ~am + carb + index)" will print dissimilarities only considering "index". Rows with NA values are omitted.
#' 
#' @examples
#' # Example with the iris dataset
#' data(iris)
#'
#' csorensendice(iris, ~Species,
#' plot = TRUE, plot_title = "Sorensen-Dice Dissimilarity Between Groups")
#'
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Example with the mtcars dataset
#' csorensendice(mtcars, ~am, plot = TRUE, 
#' plot_title = "Sorensen-Dice Dissimilarity Between Groups")
#' 
#' # Calculate the Sorensen-Dice dissimilarity for 32 car models in "mtcars" dataset
#' res <- csorensendice(mtcars, ~index)
#' 
#' @export
csorensendice <- function(dataset, formula, plot = TRUE, 
                          plot_title = "Sorensen-Dice Dissimilarity Between Groups", 
                          min_group_size = 3) {
  
  if (!is.data.frame(dataset))
    stop("The input must be a data frame.")
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0)
    stop("At least one grouping variable must be specified in the formula.")
  if (!all(grouping_vars %in% names(dataset)) && !("index" %in% grouping_vars))
    stop("Some grouping variables are not present in the dataset.")
  
  # Index mode
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In 'index' mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
    
    predictors <- setdiff(names(dataset), "index")
    
    dataset_imputed <- na.omit(dataset)
    
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    
    n <- nrow(predictors_data)
    
    D <- matrix(0, nrow = n, ncol = n)
    
    for (i in (1:n)) {
      for (j in (i:n)) {
        sim <- (2 * length(intersect(predictors_data[i,], predictors_data[j,]))) / (abs(length(predictors_data[i,])) + abs(length(predictors_data[j,])))
        
        D[i, j] <- round(1 - sim, digits = 2)
        D[j, i] <- D[i,j]
      }
    }
    rownames(D) <- colnames(D) <- as.character(dataset_imputed$index)
    
    return(list(index = list(distances = D)))
  }
  
  # Assume grouping is not "index"
  predictors <- setdiff(names(dataset), grouping_vars)
  
  if (length(predictors) != 0) {
    
    dataset_imputed <- na.omit(dataset)
    predictors_data <- model.matrix(~ . - 1, data = dataset_imputed[, predictors, drop = FALSE])
    predictors_data <- as.data.frame(predictors_data)
    variances <-  sapply(predictors_data, var)
    predictors_data <- predictors_data[, variances > 1e-9, drop = FALSE]
    predictors <- colnames(predictors_data)
    
    dataset_prepared <- cbind(dataset_imputed[grouping_vars], predictors_data)
    
    predictors_mat <- as.matrix(predictors_data)
    predictor_means <- colMeans2(predictors_mat, na.rm = TRUE)
    predictor_sds <-  apply(predictors_mat, 2, sd)
    predictor_sds[predictor_sds < 1e-9] <- 1 
    predictors_data_scaled <- scale(predictors_mat, center = predictor_means, scale = predictor_sds)
    
    dataset_prepared_scaled <- cbind(dataset_imputed[grouping_vars], as.data.frame(predictors_data_scaled))
    
    dt <- as.data.frame(dataset_prepared_scaled)
  } else if (length(predictors) == 0) {
    dataset_imputed <- na.omit(dataset)
    dataset_prepared <- (dataset_imputed)
    dt <- as.data.frame(dataset_prepared)
  }
  
  result <- list()
  plot_list <- list()
  for (grouping_var in grouping_vars) {
    groups <- split(dt, dt[[grouping_var]])
    group_sizes <-  sapply(groups, nrow)
    valid_groups <- groups[group_sizes >= min_group_size]
    
    if (length(valid_groups) < 2) {
      warning(paste("Not enough valid groups for grouping variable:", grouping_var, "- skipping."))
      next
    }
    
    group_names <- names(valid_groups)
    n <- length(valid_groups)
    
    if (length(predictors) != 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group[, predictors]))
      })
    } else if (length(predictors) == 0) {
      means <-  lapply(valid_groups, function(dt_group) {
        colMeans2(as.matrix(dt_group))
      })
    }
    
    distances <- matrix(0, nrow = n, ncol = n)
    
    rownames(distances) <- colnames(distances) <- group_names
    
    for (i in 1:n) {
      for (j in 1:n) {
        specie_1 <- means[[i]]
        specie_2 <- means[[j]]
        
        dice_similarity <- (2 * length(intersect(specie_1, specie_2))) / (abs(length(specie_1)) + abs(length(specie_2)))
        sorensendice_distance <- 1 - dice_similarity
        distances[i, j] <- round(sorensendice_distance, digits = 2)
        distances[j, i] <- distances[i, j]
      }
    }
    
    result[[grouping_var]] <- list(distances = distances)
    if (plot) {
      plot_title_i <- paste(plot_title, "-", grouping_var)
      
      dist_df <- as.matrix((distances))
      dist_df1 <- reshape2::melt(dist_df)
      
      p_obj <- ggplot2::ggplot(dist_df1, aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = plot_title, x = "", y = "", fill = "Distance") +
        ggplot2::theme_minimal()
      
      plot_list[[grouping_var]] <- p_obj
    }
  }
  
  if (plot && length(plot_list) > 0) {
    n_plots <- length(plot_list)
    n_col <- ceiling(sqrt(n_plots))
    gridExtra::grid.arrange(grobs = plot_list, ncol = n_col)
  }
  
  return(result)
}

#' @name pvaluescsore
#' @title Calculate the p_values matrix or matrices (two or more) for each factor inside variable or variables (two or more), using Sorensen-Dice dissimilarity as a base.
#' @description
#' Using the Sorensen-Dice dissimilarity for the dissimilarities calculation, this function takes a dataframe, a variable or variables (two or more), a p_value method such as "bootstrap" and "permutation" and returns the p_values matrix or matrices (two or more) between each pair of factors and a plot or plots (two or more) if the user select TRUE or leaves the parameter without argument.
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Sorensen-Dice dissimilarities matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param plot if TRUE, plot the p_values heatmap or heatmaps (two or more). The default value is TRUE.
#' @param seed Optionally, set a seed for "bootstrap" and "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A list containing a matrix or matrices (two or more) of p_values and, optionally, the plot.
#' @examples
#' # Example with the iris dataset
#' data(iris)
#' 
#' # Calculate p_values of "Species" variable in iris dataset
#' pvaluescsore(iris,~Species, pvalue.method = "bootstrap")
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Calculate p_values of "am" variable in mtcars dataset
#' pvaluescsore(mtcars,~am, 
#' pvalue.method = "permutation", seed = 134)
#' 
#' @export
pvaluescsore <- function(dataset, formula, pvalue.method = "permutation", plot = TRUE, seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "bootstrap") {
    stop("'bootstrap' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars && pvalue.method == "permutation") {
    stop("'permutation' method does not have sense with 'index'")
  }
  
  if ("index" %in% grouping_vars) {
    if (!("index" %in% names(dataset))) {
      message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
      idx <- rownames(dataset)
      is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
      if (is_default_rownames) {
        message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
        dataset$index <- seq_len(nrow(dataset))
      } else {
        message(" -> Using dataset's rownames for the index.")
        dataset$index <- idx
      }
    }
  }
  
  # Create lists to save the results (p-value matrices) and plots
  result_list <- list()
  plot_list <- list()
  excluded_message <- list()
  
  # Iterate over each grouping variable
  for (idx in seq_along(grouping_vars)) {
    grouping_var <- grouping_vars[idx]
    
    # Calculate Sorensen-Dice dissimilarities using csorensendice
    soredice_results <- csorensendice(dataset, as.formula(paste("~", grouping_var)),
                                      plot = FALSE, min_group_size = min_group_size)
    distances <- soredice_results[[grouping_var]]$distances
    
    # Number of groups
    n <- nrow(distances)
    
    # Initialize the p-values matrix
    p_values <- matrix(NA, nrow = n, ncol = n)
    rownames(p_values) <- colnames(p_values) <- rownames(distances)
    
    null_distances_list <- list() # Initialize a list to store distances
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pvalue.method == "permutation") {
      message(paste("Running permutations for", grouping_var, "..."))
      
      for (k in 1: nrow(dataset)) { # Loop for permutations
        permuted_data <- dataset
        permuted_data[[grouping_var]] <- sample(permuted_data[[grouping_var]], replace = FALSE)
        
        permuted_results <- csorensendice(permuted_data,
                                          as.formula(paste("~", grouping_var)),
                                          plot = FALSE,
                                          min_group_size = min_group_size)
        
        if (!is.null(permuted_results) && !is.null(permuted_results[[grouping_var]])) {
          null_distances_list[[k]] <- permuted_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL # Store NULL if calculation fails
        }
      }
    } else if (pvalue.method == "bootstrap") {
      message(paste("Running bootstrap samples for", grouping_var, "..."))
      for (k in 1: nrow(dataset)) { # Loop for bootstraps
        bootstrap_sample <- sample(nrow(dataset), replace = TRUE)
        bootstrap_data <- dataset[bootstrap_sample, ] 
        
        bootstrap_results <- csorensendice(bootstrap_data,
                                           as.formula(paste("~", grouping_var)),
                                           plot = FALSE,
                                           min_group_size = min_group_size)
        
        if (!is.null(bootstrap_results) && !is.null(bootstrap_results[[grouping_var]])) {
          null_distances_list[[k]] <- bootstrap_results[[grouping_var]]$distances
        } else {
          null_distances_list[[k]] <- NULL
        }
      }
    }
    
    # Now calculate p-values using the collected null distances
    for (i in 1:n) {
      for (j in i:n) {
        if (i == j) next
        
        observed_distance <- distances[i, j]
        
        # Extract the vector of distances for the cell (i,j) from all saved matrices
        replicated_distances <- sapply(null_distances_list, function(dist_matrix) {
          if (!is.null(dist_matrix) && nrow(dist_matrix) >= i && ncol(dist_matrix) >= j && !is.na(dist_matrix[i, j])) {
            return(dist_matrix[i, j])
          } else {
            return(NA)
          }
        })
        
        valid_distances <- replicated_distances[!is.na(replicated_distances)]
        
        if (length(valid_distances) == 0) {
          p_values[i, j] <- NA
        } else {
          # Calculate the p-value as a proportion of replicated distances >= to that observed
          p_values[i, j] <- round(mean(valid_distances >= observed_distance, na.rm = TRUE), digits = 2)
        }
        
        p_values[j, i] <- p_values[i, j] # Make the matrix symmetric
      }
    }
    result_list[[grouping_var]] <- p_values
    
    excluded_groups <- soredice_results[[grouping_var]]$excluded_groups
    if (length(excluded_groups) > 0 && is.null(excluded_message[[grouping_var]])) {
      excluded_message[[grouping_var]] <- paste("Grouping variable", grouping_var,
                                                "- excluded groups due to insufficient size:",
                                                paste(excluded_groups, collapse = ", "))
      message(excluded_message[[grouping_var]])
    }
    if (plot) {
      # Ensure ggplot2 is loaded and p_values_melt is created correctly
      p_values_melt <- reshape2::melt(p_values, na.rm = TRUE)
      colnames(p_values_melt) <- c("Var1", "Var2", "value")
      
      p <- ggplot2::ggplot(data = p_values_melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_colour_brewer(palette = "Greens") +
        ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
        ggplot2::geom_text(aes(label = value), color = "white", size = 4) + 
        ggplot2::labs(title = paste("p-values heatmap -", grouping_var), x = "", y = "", fill = "p_value") +
        ggplot2::theme_minimal() 
      
      plot_list[[grouping_var]] <- p
    }
  }
  if (plot && length(plot_list) > 0) {
    # Check if gridExtra is loaded
    gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  }
  return(result_list)
}

#' @name generate_report_csorensendice
#' @title Generate a Microsoft Word document about the Sorensen-Dice dissimilarity matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @description
#' This function takes a dataframe, a factor or factors (two or more) and returns a Microsoft Word document about the Sorensen-Dice dissimilarities matrix or matrices (two or more) and the p-values matrix or matrices (two or more).
#' @param dataset A dataframe.
#' @param formula A variable or variables (two or more) with factors which you want to calculate the Sorensen-Dice dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @param pvalue.method A p_value method used to calculate the matrix or matrices (two or more), the default value is "permutation". Another method is "bootstrap".
#' @param seed Optionally, set a seed for "bootstrap" or "permutation".
#' @param min_group_size Minimum group size to maintain. The default value is 3, therefore groups, inside variables, with less than 3 observations will be discarded.
#' @return A Microsoft Word document about the Sorensen-Dice dissimilarities matrix or matrices (two or more) and the p_values matrix or matrices (two or more).
#' @examples
#' # Example with iris dataset
#' data(iris)
#' 
#' # Generate a report about "Species" factor in iris dataset
#' generate_report_csorensendice(iris, ~Species, 
#' pvalue.method = 'permutation')
#' 
#' # Example with mtcars dataset
#' data(mtcars)
#' 
#' # Generate a report about "am" factor in mtcars dataset
#' generate_report_csorensendice(mtcars, ~am,
#' pvalue.method = "bootstrap", seed = 123)
#' 
#' @export
generate_report_csorensendice <- function(dataset, formula, pvalue.method = "permutation",
                                          seed = NULL, min_group_size = 3) {
  
  if (!is.data.frame(dataset)) {
    stop("The input must be a data frame.")
  }
  
  grouping_vars <- all.vars(formula)
  if (length(grouping_vars) == 0) {
    stop("At least one grouping variable must be specified in the formula.")
  }
  
  if (!all(grouping_vars %in% names(dataset)) && !identical(grouping_vars, "index")) {
    stop("Some grouping variables are not present in the dataset.")
  }
  
  if (!("index" %in% grouping_vars)) {
    # --- Index mode ---
    if ("index" %in% grouping_vars) {
      if (!("index" %in% names(dataset))) {
        message("Formula '~index' was used. In this mode, 'min_group_size' is always 1.")
        idx <- rownames(dataset)
        is_default_rownames <- is.null(idx) || all(idx == as.character(seq_len(nrow(dataset))))
        if (is_default_rownames) {
          message(" -> Using sequential numbers (1, 2, 3, ...) for the index.")
          dataset$index <- seq_len(nrow(dataset))
        } else {
          message(" -> Using dataset's rownames for the index.")
          dataset$index <- idx
        }
      }
    }  
    
    sorensendice_results <- csorensendice(dataset, formula, plot = FALSE, 
                                          plot_title = "Sorensen-Dice Dissimilarity Between Groups", 
                                          min_group_size = min_group_size)
    
    distances <- sorensendice_results
    
    
    if (pvalue.method == "permutation") {
      p_values <- pvaluescsore(dataset, formula, pvalue.method = "permutation", seed = seed, plot = FALSE, min_group_size = min_group_size)
    } else if (pvalue.method == "bootstrap") {
      p_values <- pvaluescsore(dataset, formula, pvalue.method = "bootstrap", seed = seed, plot = FALSE, min_group_size)
    }
    
    # Setup directory for figures within the temp directory (unchanged)
    dir_path <- file.path(tempdir(), "reportcsorensendice_files", "figure-docx")
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Locate the template report
    template_path <- system.file("rmarkdown", "template_report_csorensendice.Rmd", package = "cmahalanobis")
    if (file.exists(template_path)) {
      message("Template found at: ", template_path)
    } else {
      stop("Template not found! Make sure the package 'cmahalanobis' is installed and the template file exists.")
    }
    
    output_file <- file.path(tempdir(), "reportcsorensendice.docx")
    output_path <- tryCatch({
      rmarkdown::render(
        input = template_path,
        output_format = "word_document",
        output_file = output_file,
        params = list(
          distances = distances,
          p_values = p_values) # Pass the results to the template
      )
    }, error = function(e) {
      message("Error during report rendering: ", e$message)
      return(NULL)
    })
    
    if (!is.null(output_path)) {
      message("Temporary report saved in: ", output_path)
      message("If you want to store the document in a permanent path, copy it to the desired location.")
    }
    
    # Clean up temporary files
    unlink(file.path(tempdir(), "reportcsorensendice_files"), recursive = TRUE)
    temp_rmd_path <- file.path(tempdir(), basename(template_path))
    if (file.exists(temp_rmd_path)) {
      unlink(temp_rmd_path)
    }
    
    invisible(output_path)
  } else {
    stop("Report not supported with index due to permutation & bootstrap methods for p_values calculation.")
  }
}
