### Interfaces to quasi-random sequences #######################################


##' @title Korobov's sequence
##' @param n number of points (>= 2 as generator has to be in {1,..,n-1}
##' @param d dimension
##' @param generator generator in {1,..,n-1}; either a vector of length d
##'        or a single number (which is appropriately extended)
##' @param randomize character string indicating the type of randomization
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
korobov <- function(n, d = 1, generator, randomize = c("none", "shift"))
{
    stopifnot(n >= 2, d >= 1, (l <- length(generator)) == 1 || l == d,
              1 <= generator, generator <= n-1, generator %% 1 == 0)
    randomize <- match.arg(randomize)
    lim <- 2^31-1
    if(n > lim)
        stop("'n' must be <= 2^31-1")
    if(d > lim)
        stop("'d' must be <= 2^31-1")
    if(l == 1) generator <- generator^(0:(d-1)) %% n # vectorize
    u <- .Call(korobov_, n, d, generator, randomize == "shift")
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
##' @param randomize character string indicating the type of randomization
##' @param seed if provided, an integer used within set.seed() (for non-scrambling
##'        'randomize') or passed to the underlying randtoolbox::sobol() (for
##'        scrambling methods). If not provided, the non-scrambling methods respect
##'        the global .Random.seed and scrambling methods don't (they use 4711
##'        as seed)
##' @param skip the number of initial points in the sequence that are skipped
##'        (skip = 0 means the sequence starts at the origin). Note that
##'        for the scrambling methods, this simply computes n + skip points
##'        and omits the first 'skip'-many.
##' @param ... additional arguments passed to randtoolbox::sobol()
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
sobol <- function(n, d = 1, randomize = c("none", "digital.shift", "Owen",
                                          "Faure.Tezuka", "Owen.Faure.Tezuka"),
                  seed, skip = 0, ...)
{
    stopifnot(n >= 1, d >= 1, skip >= 0)
    has.seed <- hasArg(seed)
    if(is.logical(randomize))
        randomize <- if(randomize) "digital.shift" else "none" # backwards compatibility
    randomize <- match.arg(randomize)
    switch(randomize,
           "none" =, "digital.shift" = {
               if(has.seed) set.seed(seed) # to respect the seed if provided (to be compatible with scrambling methods)
               if(n > 2^31-1)
                   stop("'n' must be <= 2^31-1")
               if(d > 16510)
                   stop("'d' must be <= 16510")
               ## sobol_ <- NULL # to make CRAN check happy (for some reason not required here)
               u <- .Call(sobol_, n, d, randomize == "digital.shift", skip)
               if(d == 1) as.vector(u) else u
           },
           "Owen" =, "Faure.Tezuka" =, "Owen.Faure.Tezuka" = {
               scrambling <- if(randomize == "Owen") 1 else if (randomize == "Faure.Tezuka") 2 else 3
               if(!has.seed) seed <- 4711 # randtoolbox::sobol's default
               res <- randtoolbox::sobol(n + skip, dim = d, scrambling = scrambling,
                                         seed = seed, ...)
               if(d == 1) res[(1+skip):(n+skip)] else res[(1+skip):(n+skip),]
           }, stop("Wrong 'randomize'"))
}
