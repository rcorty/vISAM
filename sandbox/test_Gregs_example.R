library(wISAM)

d <- readRDS('sandbox/test_inputs_from_Greg.RDS')
str(d)

g <- GenomeScan$new(y = d$y, X = d$X, K = d$K)

g <- g$prep_scan()

g$Intermediates_per_scan$LL_null
g$Intermediates_per_scan$h2_null

g$Results$h2


library(emma)

e <- emma.ML.LRT(ys = matrix(data = d$y, ncol = 1),
                 xs = matrix(data = c(0, rep(1, 280)), nrow = 281, ncol = 1),
                 K = d$K)

e$ML1s
e$ML0s
e$vgs/(e$vgs + e$ves)
