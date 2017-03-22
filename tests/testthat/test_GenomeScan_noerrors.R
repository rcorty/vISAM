context('Testing whether GenomeScan runs without error')

library(MASS)

set.seed(seed = 27599)

n <- 1e2
c <- 4
p <- 1e2
s <- sample(x = c(-2, -1, 0, 1, 2), size = n, replace = TRUE) + 10
w <- 3/sqrt(s)

X <- matrix(data = sample(x = n*c), nrow = n)
G <- matrix(data = rbinom(n = n*p, size = 2, prob = 0.3), nrow = n, ncol = p)
K <- 1 - as.matrix(x = dist(x = G, method = 'manhattan'))/p

A <- mvrnorm(n = 1,
             mu = rep(0, n),
             # mu = G %*% replace(x = rep(0, p), list = 1, values = 1),
             Sigma = K)
E <- rnorm(n = n, sd = w)
y <- A + E

gs <- GenomeScan$new(y = y, X = X, G = G, K = K)
wgs <- GenomeScan$new(y = y, X = X, G = G, K = K, w = w)



test_that(desc = 'Create a GenomeScan object',
          code = {

              # valid inputs
              expect_is(object = GenomeScan$new(y = y, X = X, G = G, K = K),
                        class = 'GenomeScan')
              expect_is(object = GenomeScan$new(y = y, X = X, G = G, K = K, w = w),
                        class = 'GenomeScan')

              # one or more necessary inputs missing
              expect_error(object = GenomeScan$new())
              expect_error(object = GenomeScan$new(       X = X, G = G, K = K))
              expect_error(object = GenomeScan$new(y = y,        G = G, K = K))
              expect_error(object = GenomeScan$new(y = y, X = X,        K = K))

              # right inputs, one with the wrong dimension
              expect_error(object = GenomeScan$new(y = y[-1], X = X,      G = G,      K = K,      w = w))
              expect_error(object = GenomeScan$new(y = y,     X = X[-1,], G = G,      K = K,      w = w))
              expect_error(object = GenomeScan$new(y = y,     X = X,      G = G[-1,], K = K,      w = w))
              expect_error(object = GenomeScan$new(y = y,     X = X,      G = G,      K = K[-1,], w = w))
              expect_error(object = GenomeScan$new(y = y,     X = X,      G = G,      K = K,      w = w[-1]))

          })


test_that(desc = 'Prepare for a scan',
          code = {

              expect_is(object = gs$prep_scan(),  class = 'GenomeScan')
              expect_is(object = wgs$prep_scan(), class = 'GenomeScan')

          })


test_that(desc = 'Conduct a scan',
          code = {

              expect_is(object = gs$conduct_scan(),  class = 'GenomeScan')
              expect_is(object = wgs$conduct_scan(), class = 'GenomeScan')

          })