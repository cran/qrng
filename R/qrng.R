### Interfaces to quasi-random sequences #######################################


##' @title Korobov's sequence
##' @param n number of points (>= 2 as generator has to be in {1,..,n-1}
##' @param d dimension
##' @param generator generator in {1,..,n-1}; either a vector of length d
##'        or a single number (which is appropriately extended)
##' @param randomize logical indicating whether the point set should be
##'        randomized
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
korobov <- function(n, d = 1, generator, randomize = FALSE)
{
    stopifnot(n >= 2, d >= 1, (l <- length(generator)) == 1 || l == d,
              1 <= generator, generator <= n-1, generator %% 1 == 0)
    lim <- 2^31-1
    if(n > lim)
        stop("'n' must be <= 2^31-1")
    if(d > lim)
        stop("'d' must be <= 2^31-1")
    if(l == 1) generator <- generator^(0:(d-1)) %% n # vectorize
    u <- .Call(korobov_, n, d, generator, randomize)
    if(d == 1) as.vector(u) else u
}

##' @title Generalized Halton sequence
##' @param n number of points
##' @param d dimension
##' @param method character string indicating which sequence is generated
##'        (generalized Halton or (plain) Halton)
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
ghalton <- function(n, d = 1, method = c("generalized", "halton"))
{
    stopifnot(n >= 1, d >= 1)
    method <- match.arg(method)
    if(n > 2^32-1)
        stop("'n' must be <= 2^32-1")
    if(d > 360)
        stop("'d' must be <= 360")
    ## ghalton_ <- NULL # to make CRAN check happy (for some reason not required here)
    u <- .Call(ghalton_, n, d, method)
    if(d == 1) as.vector(u) else u
}

##' @title Sobol sequence
##' @param n number of points
##' @param d dimension
##' @param randomize logical indicating whether a digital shift should be
##'        included
##' @param skip number of initial terms in the sequence that are skipped
##'        (skip = 0 means the sequence starts with the origin)
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
sobol <- function(n, d = 1, randomize = FALSE, skip = 0)
{
    stopifnot(n >= 1, d >= 1, is.logical(randomize), skip >= 0)
    if(n > 2^31-1)
        stop("'n' must be <= 2^31-1")
    if(d > 16510)
        stop("'d' must be <= 16510")
    ## sobol_ <- NULL # to make CRAN check happy (for some reason not required here)
    u <- .Call(sobol_, n, d, randomize, skip)
    if(d == 1) as.vector(u) else u
}

