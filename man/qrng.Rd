\name{qrng}
\alias{korobov}
\alias{ghalton}
\alias{sobol}
\title{Compute Quasi-Random Sequences}
\description{
  Computing Korobov, generalize Halton and Sobol quasi-random sequences.
}
\usage{
korobov(n, d, generator, randomize = FALSE)
ghalton(n, d, method = c("generalized", "halton"))
sobol(n, d, randomize = FALSE)
}
\arguments{
  \item{n}{Number \eqn{n} of points to be generated \eqn{\ge 2}{>= 2}.}
  \item{d}{Dimension \eqn{d}.}
  \item{generator}{A \code{\link{numeric}} of length \eqn{d} or length 1
    (in which case it is appropriately extended to length \eqn{d}). All numbers
    must be in \eqn{\{1,\dots,n\}}{{1,...,n}} and must be (coercible to)
    integers.}
  \item{randomize}{A \code{\link{logical}} indicating whether the point
    set should be randomized (for \code{sobol()}: a digital shift).}
  \item{method}{A \code{\link{character}} string indicating which
    sequence is generated, generalized Halton or (plain) Halton.}
}
\value{
  \code{korobov()} and \code{ghalton()} return an
  \eqn{(n,d)}-\code{\link{matrix}}; for \eqn{d=1} an \eqn{n}-vector
  is returned.
}
\details{
  Note that these procedures call fast C code. The following
  restrictions apply:
  \describe{
    \item{korobov()}{\code{n},\code{d} must be \eqn{\le 2^{31}-1}{<= 2^31-1}.}
    \item{ghalton()}{\code{n} must be \eqn{\le 2^{32}-1}{<= 2^32-1} and
      \code{d} must be \eqn{\le 360}{<= 360}.}
    \item{sobol()}{\code{n} must be \eqn{\le 2^{31}-1}{<= 2^31-1} and
      \code{d} must be \eqn{\le 360}{<= 360}.}
  }

  The choice of parameters for \code{korobov()} is crucial for the quality of
  this quasi-random sequence (only basic sanity checks are
  conducted). For more details, see l'Ecuyer and Lemieux (2000).

  The generalized Halton sequence uses the scrambling factors of
  Faure and Lemieux (2009).
}
\author{Marius Hofert and Christiane Lemieux}
\references{
  Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
  A Comparative Study. \emph{ACM-TOMACS} \bold{19}(4), Article 15.

  l'Ecuyer, P., Lemieux, C. (2000). Variance Reduction via Lattice
  Rules. \emph{Stochastic Models and Simulation}, 1214--1235.

  Lemieux, C., Cieslak, M., Luttmer, K. (2004). RandQMC User's guide.
  See https://www.math.uwaterloo.ca/~clemieux/randqmc/guide.pdf
}
%% \seealso{
%% }
\examples{
n <- 1021 # prime
d <- 4 # dimension

## Korobov's sequence
generator <- 76 # see l'Ecuyer and Lemieux
u <- korobov(n, d = d, generator = generator)
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Randomized Korobov's sequence
set.seed(271)
u <- korobov(n, d = d, generator = generator, randomize = TRUE)
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Generalized Halton sequence (randomized by definition)
set.seed(271)
u <- ghalton(n, d)
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Randomized Sobol sequence (with digital shift)
set.seed(271)
u <- sobol(n, d, randomize = TRUE)
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))
}
\keyword{distribution}