\name{test_functions}
\alias{sum_of_squares}
\alias{sobol_g}
\alias{exceedance}
\title{Test Functions}
\description{
  Functions for testing low-discrepancy sequences.
}
\usage{
sum_of_squares(u)
sobol_g(u, copula = indepCopula(dim = ncol(u)), alpha = 1:ncol(u), ...)
exceedance(x, q, p = 0.99, method = c("indicator", "individual.given.sum.exceeds",
                                      "sum.given.sum.exceeds"))
}
\arguments{
  \item{u}{\eqn{(n, d)}-matrix containing \eqn{n}
    \eqn{d}-dimensional realizations (of a potential quasi-random number
    generator). For \code{sum_of_squares()} these need to be marginally
    standard uniform and for \code{sobol_g()} they need to follow the
    copula specified by \code{copula}.}
  \item{copula}{\code{\link{Copula}} object for which the inverse
    Rosenblatt transformation exists.}
  \item{alpha}{vector of parameters of Sobol's g test function.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{cCopula}()}.}
  \item{x}{\eqn{(n, d)}-matrix containing \eqn{n} \eqn{d}-dimensional
    realizations.}
  \item{q}{
    \describe{
      \item{"indicator"}{\eqn{d}-vector containing the componentwise
	thresholds; if a number it is recycled to a \eqn{d}-vector.}
      \item{"individual.given.sum.exceeds", "sum.given.sum.exceeds"}{
        threshold for the sum (row sums of \code{x}).}
    }
  }
  \item{p}{
    If \code{q} is not provided, the probability \code{p} is used to
    determine \code{q}.
    \describe{
      \item{"indicator"}{\eqn{d}-vector containing the probabilities
	determining componentwise thresholds via empirical quantiles;
        if a number, it is recycled to a \eqn{d}-vector.}
      \item{"individual.given.sum.exceeds", "sum.given.sum.exceeds"}{
	probability determining the threshold for the sum (row sums of
	\code{x}) via the corresponding empirical quantile.}
    }
  }
  \item{method}{\code{\link{character}} string indicating the type of
    exceedance computed (see Section Value below).}
}
\value{
  \code{sum_of_squares()} returns an \eqn{n}-vector
  (\code{\link{numeric}(n)}) with the rowwise computed scaled sum
  of squares (theoretically integrating to 1).

  \code{sobol_g()} returns an \eqn{n}-vector (\code{\link{numeric}(n)})
  with the rowwise computed Sobol' g functions.

  \code{exceedance()}'s return value depends on \code{method}:
  \describe{
    \item{"indicator"}{returns indicators whether,
      componentwise, \code{x} exceeds the threshold determined by \code{q}.}
    \item{"individual.given.sum.exceeds"}{returns all rows of \code{x}
      whose sum exceeds the threshold determined by \code{q}.}
    \item{"sum.given.sum.exceeds"}{returns the row sums of those
      rows of \code{x} whose sum exceeds the threshold determined by \code{q}.}
  }
}
\details{
  For examples see the demo \code{man_test_functions}.

  See \code{ES_np(<matrix>)} from \pkg{qrmtools} for another test function.
}
\author{Marius Hofert and Christiane Lemieux}
\references{
  Radovic, I., Sobol', I. M. and Tichy, R. F. (1996). Quasi-Monte Carlo
  methods for numerical integration: Comparison of different low
  discrepancy sequences. \emph{Monte Carlo Methods and Applications}
  \bold{2}(1), 1--14.

  Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
  A Comparative Study. \emph{ACM-TOMACS} \bold{19}(4), Article 15.

  Owen, A. B. (2003). The dimension distribution and quadrature test
  functions. \emph{Stat. Sinica} \bold{13}, 1-–17.

  Sobol', I. M. and Asotsky, D. I. (2003). One more experiment on
  estimating high-dimensional integrals by quasi-Monte Carlo methods.
  \emph{Math. Comput. Simul.} \bold{62}, 255–-263.
}
\examples{
## Generate some (here: copula, pseudo-random) data
library(copula)
set.seed(271)
cop <- claytonCopula(iTau(claytonCopula(), tau = 0.5)) # Clayton copula
U <- rCopula(1000, copula = cop)

## Compute sum of squares test function
mean(sum_of_squares(U)) # estimate of E(3(sum_{j=1}^d U_j^2)/d)

## Compute the Sobol' g test function
if(packageVersion("copula") >= "0.999-20")
    mean(sobol_g(U)) # estimate of E(<Sobol's g function>)

## Compute an exceedance probability
X <- qnorm(U)
mean(exceedance(X, q = qnorm(0.99))) # fixed threshold q
mean(exceedance(X, p = 0.99)) # empirically estimated marginal p-quantiles as thresholds

## Compute 99% expected shortfall for the sum
mean(exceedance(X, p = 0.99, method = "sum.given.sum.exceeds"))
## Or use ES_np(X, level = 0.99) from 'qrmtools'
}
\keyword{utilities}