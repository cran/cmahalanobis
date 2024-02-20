# testthat.R

# Load the testthat package
library(testthat)

# Load your package
library(cmahalanobis)

# Run the tests in the tests/testthat folder
test_check("cmahalanobis")
