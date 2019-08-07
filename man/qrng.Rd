\name{qrng}
\alias{korobov}
\alias{ghalton}
\alias{sobol}
\title{Compute Quasi-Random Sequences}
\description{
  Computing Korobov, generalize Halton and Sobol' quasi-random sequences.
}
\usage{
korobov(n, d = 1, generator, randomize = c("none", "shift"))
ghalton(n, d = 1, method = c("generalized", "halton"))
sobol  (n, d = 1, randomize = c("none", "digital.shift", "Owen", "Faure.Tezuka",
                                "Owen.Faure.Tezuka"), seed, skip = 0, ...)
}
\arguments{
  \item{n}{Number \eqn{n} of points to be generated \eqn{\ge 2}{>= 2}.}
  \item{d}{Dimension \eqn{d}.}
  \item{generator}{A \code{\link{numeric}} of length \eqn{d} or length 1
    (in which case it is appropriately extended to length \eqn{d}). All numbers
    must be in \eqn{\{1,\dots,n\}}{{1,...,n}} and must be (coercible to)
    integers.}
  \item{randomize}{A \code{\link{character}} string indicating the type
    of randomization for the point set
    \describe{
      \item{"none"}{no randomization.}
      \item{"shift"}{uniform random variate modulo 1.}
      \item{"digital.shift"}{a digital shift.}
      \item{"Owen", "Faure.Tezuka","Owen.Faure.Tezuka"}{calls
	\code{\link{sobol}()} from package \pkg{randtoolbox} with
        \code{scrambling} being 1, 2 and 3, respectively.}
      }
      If \code{randomize} is a \code{\link{logical}}, then
      it is interpreted as \code{"none"} if \code{FALSE}
      and \code{"digital.shift"} if \code{TRUE}.
  }
  \item{method}{A \code{\link{character}} string indicating which
    sequence is generated, generalized Halton or (plain) Halton.}
  \item{seed}{if provided, an integer used within
    \code{\link{set.seed}()} for the non-scrambling \code{randomize}
    methods (so \code{"none"} or \code{"digital.shift"}) or passed
    to the underlying \code{\link{sobol}()} from package
    \pkg{randtoolbox} for the scrambling methods. If not provided,
    the non-scrambling methods respect a global \code{\link{set.seed}()}
    but scrambling methods do not (they then use 4711 as default);
    see \code{\link{sobol}()} from package \pkg{randtoolbox}.}
  \item{skip}{number of initial points in the sequence to be skipped
    (\code{skip = 0} means the sequence starts with at the origin).
    Note that for the scrambling methods this simply computes \code{n +
      skip} points and omits the first \code{skip}-many.}
  \item{\dots}{additional arguments passed to \code{\link{sobol}()}
    from package \pkg{randtoolbox}.}
}
\value{
  \code{korobov()} and \code{ghalton()} return an
  \eqn{(n,d)}-\code{\link{matrix}}; for \eqn{d=1} an \eqn{n}-vector
  is returned.
}
\details{
  For \code{sobol()} examples see \code{demo(sobol_examples)}.
  In particular, be careful when using
  \code{skip > 0} and \code{randomize = TRUE}; in this case, choosing
  a wrong seed (or no seed) might lead to a bad sequence.

  Note that these procedures call fast C code. The following
  restrictions apply:
  \describe{
    \item{korobov()}{\code{n},\code{d} must be \eqn{\le 2^{31}-1}{<= 2^31-1}.}
    \item{ghalton()}{\code{n} must be \eqn{\le 2^{32}-1}{<= 2^32-1} and
      \code{d} must be \eqn{\le 360}{<= 360}.}
    \item{sobol()}{if \code{randomize = "none"} or \code{randomize =
      "digital.shift"}, \code{n} must be \eqn{\le 2^{31}-1}{<= 2^31-1} and
      \code{d} must be \eqn{\le 16510}{<= 16510}.}
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
u <- korobov(n, d = d, generator = generator, randomize = "shift")
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Generalized Halton sequence (randomized by definition)
set.seed(271)
u <- ghalton(n, d)
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## For sobol() examples, see demo(sobol_examples)
}
\keyword{distribution}