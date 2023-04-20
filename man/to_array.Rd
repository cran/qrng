\name{to_array}
\alias{to_array}
\title{Compute Matrices to Arrays}
\description{
  Converting higher-dimensional matrices of quasi-random numbers to
  arrays of specific formats.
}
\usage{
to_array(x, f, format = c("(n*f,d)", "(n,f,d)"))
}
\arguments{
  \item{x}{(\eqn{n},\eqn{f}\eqn{d})-matrix of quasi-random numbers to
    be converted.}
  \item{f}{factor \eqn{f\ge 1}{f >= 1} dividing \code{ncol{x}}.}
  \item{format}{\code{\link{character}} string indicating the output format
    to which \code{x} should be converted.}
}
\value{
  \code{(n * f, d)-\link{matrix}} or \code{(n, f, d)-\link{array}}
  depending on the chosen \code{format}.
}
\details{
  \code{to_array()} is helpful for converting quasi-random numbers
  to time series paths.
}
\author{Marius Hofert}
\seealso{
  \code{\link{korobov}()}, \code{\link{ghalton}()}, \code{\link{sobol}()}.
}
\examples{
## Basic call
N <- 4 # replications
n <- 3 # time steps
d <- 2 # dimension
set.seed(271) # note: respected for the choice of 'randomize'
x <- sobol(N, d = n * d, randomize = "digital.shift") # higher-dim. Sobol'
stopifnot(dim(to_array(x, f = n)) == c(N * n, d)) # conversion and check
stopifnot(dim(to_array(x, f = n, format = "(n,f,d)")) == c(N, n, d))

## See how the conversion is done
(x <- matrix(1:(N * n * d), nrow = N, byrow = TRUE))
to_array(x, f = n) # => (n * d)-column x was blocked in n groups of size d each
}
\keyword{utilities}