library(SimHaploidPop)
library(lattice)

emma_eigen <- function(Z, K) {
    eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE)
    return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors))))
}

reconstitute <- function(e) {
    return(e$vectors %*% diag(e$values) %*% t(e$vectors))
}

is_I <- function(m) {
    sum(abs(m - diag(x = nrow(m)))) < 1e-6
}

set.seed(1)

num.strains <- 20
num.indivs.per.strain <- sample(x = c(1, 2, 3), #c(5, 10, 15),
                                size = num.strains, replace = TRUE)

D <- diag(x = rep(1, num.indivs))

num.indivs <- sum(num.indivs.per.strain)
strain.by.indiv <- rep(1:num.strains, times = num.indivs.per.strain)
Z <- model.matrix(~ 0 + factor(strain.by.indiv))

d <- dist(SimHaploidPop(stop.at.pop.size = num.strains),
          method = 'manhattan',
          diag = TRUE, upper = TRUE) %>% as.matrix()


k <- (max(d) - d)/(max(d) - min(d))
K <- Z %*% k %*% t(Z)

system.time(eigen_k <- eigen(k))
system.time(eigen_K <- eigen(K))
system.time(lr_eigen_K <- emma_eigen(Z = Z, K = k))

str(eigen_k)
str(eigen_K)
str(lr_eigen_K)

system.time(message(mean(reconstitute(e = eigen_k) - k)))
system.time(message(mean(reconstitute(e = eigen_K) - K)))
system.time(message(mean(reconstitute(e = lr_eigen_K) - K)))

is_I(m = t(eigen_k$vectors) %*% eigen_k$vectors)
is_I(m = (eigen_k$vectors) %*% t(eigen_k$vectors))

is_I(m = t(eigen_K$vectors) %*% eigen_K$vectors)
is_I(m = (eigen_K$vectors) %*% t(eigen_K$vectors))

is_I(m = t(lr_eigen_K$vectors) %*% (lr_eigen_K$vectors))
is_I(m = (lr_eigen_K$vectors) %*% t(lr_eigen_K$vectors))


lr_eigen_K_aug <- list(values = c(lr_eigen_K$values, rep(0, num.indivs - num.strains)),
                       vectors = cbind(lr_eigen_K$vectors, matrix(data = 0, nrow = num.indivs, ncol = num.indivs - num.strains)))

mean(reconstitute(e = lr_eigen_K_aug) - K)

is_I(m = t(lr_eigen_K_aug$vectors) %*% (lr_eigen_K_aug$vectors))
is_I(m = (lr_eigen_K_aug$vectors) %*% t(lr_eigen_K_aug$vectors))






is_I(m = (D %*% lr_eigen_K$vectors) %*% t(D %*% lr_eigen_K$vectors))
is_I(m = t(D %*% lr_eigen_K$vectors) %*% (D %*% lr_eigen_K$vectors))



a <- diag(1/sqrt(diag(D))) %*% Z %*% eigen_k$vectors
is_I(a %*% t(a))
is_I(t(a) %*% (a))

b <- t(eigen_k$vectors) %*% t(Z) %*% diag(1/sqrt(diag(D)))
is_I(b %*% t(b))
is_I(t(b) %*% b)

