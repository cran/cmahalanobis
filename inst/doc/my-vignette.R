## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 4,
  fig.height = 4,
  out.width = '40%'
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages(”cmahalanobis”)

## -----------------------------------------------------------------------------
library(cmahalanobis)
num_observations <- 100
num_variables <- 5
group1 <- data.frame(a = c(1,2,NA,4), b = c(3,4,5,NA))
group2 <- data.frame(a = c(2,5,NA,3), b = c(3,4,5,5))
groups <- list(group1, group2)
distances <- cmahalanobis(groups, plot = TRUE, p.value = TRUE)
distances

## -----------------------------------------------------------------------------
group1 <- data.frame(a = c(1,2,NaN,4), b = c(3,4,5,NaN))
group2 <- data.frame(a = c(2,5,NaN,3), b = c(3,4,5,5))
groups <- list(group1, group2)
distances <- cmahalanobis(groups, plot = FALSE, p.value = TRUE)
distances

## -----------------------------------------------------------------------------
iris_list <- split(iris, iris$Species)

## -----------------------------------------------------------------------------
# Print iris_list
iris_list

## -----------------------------------------------------------------------------
res <- cmahalanobis(iris_list, p.value = TRUE)

## -----------------------------------------------------------------------------
res

## -----------------------------------------------------------------------------
# Create a dataframe where only ”am = 0” is present
auto <- subset(mtcars, am == 0)
# Remove the variable ”am = 0”
auto <- auto [, -9]
# Create a dataframe where only ”am = 1” is present
manual <- subset(mtcars, am == 1)
# Remove the variable ”am = 1”
manual <- manual[, -9]
# Create a list with the two groups of cars
groups <- list(auto, manual)

## -----------------------------------------------------------------------------
# Print groups
groups

## -----------------------------------------------------------------------------
res <- cmahalanobis(groups, plot = TRUE, p.value = TRUE)
res

## -----------------------------------------------------------------------------
# Load cmahalanobis package
library(cmahalanobis)
# Define the number of observations and variables for each groups
num_observations <- 100
num_variables <- 5
# We generate three groups of simulated data with normal distribution
set.seed(123) # For the reproducibility of results
group1 <- as.data.frame(matrix(rnorm(num_observations * num_variables), 
                               nrow = num_observations))
group2 <- as.data.frame(matrix(rnorm(num_observations * num_variables), 
                               nrow = num_observations))
group3 <- as.data.frame(matrix(rnorm(num_observations * num_variables), 
                               nrow = num_observations))
# Create a list of three groups of data
groups <- list(group1, group2, group3)
# Calculate Mahalanobis distance with cmahalanobis function
distances <- cmahalanobis(groups, p.value = TRUE)
# Visualize the distance matrix
distances

