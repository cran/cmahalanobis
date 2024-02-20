% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmahalanobis.R
\name{cmahalanobis}
\alias{cmahalanobis}
\title{Calculate the Mahalanobis distance for each species}
\usage{
cmahalanobis(dataset)
}
\arguments{
\item{dataset}{A list of data frames, where each data frame contains the observations of a species.}
}
\value{
A matrix with the Mahalanobis distances between each pair of species.
}
\description{
. This function takes a list of data frames as input, where each data frame contains the observations of a species, and returns a matrix with the Mahalanobis distances between each pair of species.
}
\examples{
# Example with the iris dataset
library(stats)
# Split the data into 3 parts for each species
setosa <- subset(iris, Species == "setosa")
setosa <- setosa[,-5]
versicolor <- subset(iris, Species == "versicolor")
versicolor <- versicolor[,-5]
virginica <- subset(iris, Species == "virginica")
virginica <- virginica[,-5]

# Create a list with the three groups of flowers
groups <- list(setosa, versicolor, virginica)

# Calculate the Mahalanobis distance with the cmahalanobis function
cmahalanobis(groups)

# Example with the mtcars dataset
library(stats)
# Split the data into 2 parts for each type of transmission
auto <- subset(mtcars, am == 0)
auto <- auto[,-9]
manual <- subset(mtcars, am == 1)
manual <- manual[,-9]

# Create a list with the two groups of cars
groups <- list(auto, manual)

# Calculate the Mahalanobis distance with the cmahalanobis function
cmahalanobis(groups)

}