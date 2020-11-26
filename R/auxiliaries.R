### Auxiliaries ################################################################

##' @title Convert Matrices to Arrays
##' @param x (n, f * d)-matrix, e.g. sobol(n, d = f * d, randomize = "digital.shift")
##' @param f integer f >= 1 dividing ncol(x)
##' @param format string specifying how output should look like
##' @return an array of the form as given by 'format'
##' @author Marius Hofert
##' @note This is helpful for time series applications where f = number of time
##'       steps to simulate
to_array <- function(x, f, format = c("(n*f,d)", "(n,f,d)"))
{
    if(!is.matrix(x))
        x <- rbind(x) # convert vectors to 1-row matrices
    stopifnot(f %% 1 == 0, f >= 1)
    if(f == 1) return(x)
    dm <- dim(x)
    n  <- dm[1]
    d. <- dm[2] # = f * d
    if(d. %% f != 0)
        stop("'f' must divide d = ncol(x)")
    d <- d./f
    format <- match.arg(format)
    switch(format,
           "(n*f,d)" = { # convert (n, f * d)-matrix into (n * f, d)-matrix
               matrix(t(x), ncol = d, byrow = TRUE)
           },
           "(n,f,d)" = { # convert (n, f * d)-matrix into (n, f, d)-array
               array(x, dim = c(n, f, d))
           },
           stop("Wrong 'format'"))
}
