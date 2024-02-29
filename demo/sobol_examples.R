## By Marius Hofert

## Examples of sobol() calls (presented here instead of on ?sobol) to avoid
## a non-reproducible triggering of valgrind (2019-08-06).

library(qrng)

n <- 1021 # prime
d <- 4 # dimension

## Randomized Sobol' sequence (with digital shift)
set.seed(271)
u <- sobol(n, d, randomize = "digital.shift")
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Randomized Sobol' sequence (with Owen scrambling)
## Note: This is calling spacefillr::generate_sobol_owen_set()
u <- sobol(n, d, randomize = "Owen") # auto-generates a seed
pairs(u, gap = 0, pch = ".", labels = as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

## Check whether a Sobol' sequence of size 2*n equals one of size n
## and, concatenated, one of size n with the first n numbers skipped
f <- function(n)
{
    set.seed(271)
    a <- sobol(2*n, randomize = "digital.shift")
    set.seed(271)
    b1 <- sobol(n, randomize = "digital.shift")
    set.seed(271)
    b2 <- sobol(n, randomize = "digital.shift", skip = n)
    all(all.equal(apply(cbind(a, c(b1, b2)), 1, diff), rep(0, 2*n)))
}
stopifnot(sapply(1:10, f)) # check for n = 1, ..., 10

## Careful when using skip > 0 and randomization => seed matters!

## Drawing all points at once (works, of course)
n <- 32
set.seed(5)
U.2n <- sobol(2*n, d = 2, randomize = "digital.shift")
plot(U.2n, main = "All points generated at once",
     xlab = expression(U[1]), ylab = expression(U[2]))

## Drawing successively with the same seed (typically the right approach)
set.seed(5)
U.n.1 <- sobol(n, d = 2, randomize = "digital.shift")
set.seed(5) # => same seed
U.n.2 <- sobol(n, d = 2, randomize = "digital.shift", skip = n)
U.2n.same.seed <- rbind(U.n.1, U.n.2)
plot(U.2n.same.seed,
     main = "Drawing successively, with the same seed",
     xlab = expression(U[1]), ylab = expression(U[2]))
stopifnot(all.equal(U.2n, U.2n.same.seed)) # sanity check

## Drawing successively but with different seeds (typically the wrong approach)
set.seed(5)
U.n.1 <- sobol(n, d = 2, randomize = "digital.shift", skip = 0)
set.seed(22)
U.n.2 <- sobol(n, d = 2, randomize = "digital.shift", skip = n)
U.2n.different.seed <- rbind(U.n.1, U.n.2)
plot(U.2n.different.seed,
     main = "Drawing successively, with different seeds",
     xlab = expression(U[1]), ylab = expression(U[2]))
