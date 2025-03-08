% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmahalanobis.R
\name{cbraycurtis}
\alias{cbraycurtis}
\title{Calculate the Bray-Curtis distance for each species.}
\usage{
cbraycurtis(
  dataset,
  formula,
  plot = TRUE,
  plot_title = "Bray-Curtis Distance Between Groups"
)
}
\arguments{
\item{dataset}{A dataframe.}

\item{formula}{A factor which you want to calculate the Bray-Curtis distances matrix.}

\item{plot}{Logical, if TRUE, a plot of Bray-Curtis distances matrix is displayed.}

\item{plot_title}{The title to be used for the plot if plot is TRUE.}
}
\value{
A matrix containing Bray-Curtis distances between each pair of groups and the plot.
}
\description{
This function takes a dataframe and a factor in input, and returns a matrix with the Bray-Curtis distances about it.
}
\examples{
# Example with the iris dataset

cbraycurtis(iris, ~Species, plot = TRUE, plot_title = "Bray-Curtis Distance Between Groups")

# Example with the mtcars dataset
cbraycurtis(mtcars, ~am, plot = TRUE, plot_title = "Bray-Curtis Distance Between Groups")

}
