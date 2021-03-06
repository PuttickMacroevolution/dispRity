% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho.utilities.R
\name{get.contrast.matrix}
\alias{get.contrast.matrix}
\title{Generates a contrast matrix.}
\usage{
get.contrast.matrix(matrix)
}
\arguments{
\item{matrix}{a discrete morphological character matrix.}
}
\description{
Creates a contrast matrix using the observed character states in an input matrix.
}
\examples{
## A random multistate matrix
matrix <- matrix(sample(c(0,1,2), 100, TRUE), 10, 10)

## Get the contrast matrix
get.contrast.matrix(matrix)

## Adding inapplicable and missing data to the matrix
matrix[sample(1:100, 10)] <- "?"
matrix[sample(1:100, 10)] <- "-"

## Get the contrast matrix
get.contrast.matrix(matrix)

}
\seealso{
\code{\link{check.morpho}}
}
\author{
Thomas Guillerme
}
