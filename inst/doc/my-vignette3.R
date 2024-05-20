## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(cmahalanobis)
# Load iris dataset
data(iris)
# Split data into three parts
setosa <- subset(iris, Species == "setosa")
setosa <- setosa[,-5] # Remove the column of specie
versicolor <- subset(iris, Species == "versicolor")
versicolor <- versicolor[,-5] # Remove the column of specie
virginica <- subset(iris, Species == "virginica")
virginica <- virginica[,-5] # Remove the column of specie
# Create a list with the three groups of flowers
groups <- list(setosa, versicolor, virginica)
cmahalanobis(groups, plot = TRUE, p.value = TRUE)

## -----------------------------------------------------------------------------
ceuclide(groups, plot = TRUE, p.value = TRUE, num.permutations = 10)

## -----------------------------------------------------------------------------
cmanhattan(groups, plot = TRUE, p.value = TRUE, num.permutations = 10)

## -----------------------------------------------------------------------------
cchebyshev(groups, plot = TRUE, p.value = TRUE, num.permutations = 50)

## -----------------------------------------------------------------------------
# Split the data into 2 parts for each type of transmission
auto <- subset(mtcars, am == 0)
auto <- auto[,-9]
manual <- subset(mtcars, am == 1)
manual <- manual[,-9]

# Create a list with the two groups of cars
groups <- list(auto, manual)
cmahalanobis(groups, plot = TRUE, p.value = TRUE)

## -----------------------------------------------------------------------------
ceuclide(groups, plot = TRUE, p.value = TRUE, num.permutations = 10)

## -----------------------------------------------------------------------------
cmanhattan(groups, plot = TRUE, p.value = TRUE, num.permutations = 10)

## -----------------------------------------------------------------------------
cchebyshev(groups, plot = TRUE, p.value = TRUE, num.permutations = 50)

## -----------------------------------------------------------------------------
# Load cmahalanobis package
library(cmahalanobis)
# Define the number of observations and variables for each groups
num_observations <- 100
num_variables <- 5
# We generate three groups of simulated data with normal distribution
set.seed(123) # For the reproducibility of results
group1 <- as.data.frame(matrix(rnorm(num_observations * num_variables), nrow = num_observations))
group2 <- as.data.frame(matrix(rnorm(num_observations * num_variables), nrow = num_observations))
group3 <- as.data.frame(matrix(rnorm(num_observations * num_variables), nrow = num_observations))
# Create a list of three groups of data
groups <- list(group1, group2, group3)
# Calculate Mahalanobis distance with cmahalanobis function
cmahalanobis(groups, plot = TRUE, p.value = TRUE)

## -----------------------------------------------------------------------------
ceuclide(groups, plot = TRUE, p.value = TRUE, num.permutations = 190)

## -----------------------------------------------------------------------------
cmanhattan(groups, plot = TRUE, p.value = TRUE, num.permutations = 10)

## -----------------------------------------------------------------------------
cchebyshev(groups, plot = TRUE, p.value = TRUE, num.permutations = 50)

