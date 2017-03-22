GenomeScan <-
    setRefClass(
        Class = 'GenomeScan',
        fields = c('Data',                    # y, X, G, K, and weights -- things that have interpretable meaning
                   'Intermediates_per_scan',  # L, eigen_L, LL_0 -- uninterpretable computational tools that are constant for a GenomeScan
                   'Intermediates_per_locus', # M -- uninterpretable computational tools that are constant for a locus
                   'Intermediates_per_fit',   # sigma_a and sigma_e -- interpretable, but changing very fast, don't look here for results
                   'Results'),                # maximum likelihood LL, sigma_a, and sigma_e for each locus

        methods = list(

            #   INPUTS
            #   y - vector of length n - the phenotype of each of n individuals
            #   X - matrix of dimension n-by-c - the covariate value of each individual for c covariates
            #   G - matrix of dimension n-by-p - the genotype of each individual at p loci
            #   K - matrix of dimension n-by-n - the covariance of the phenotype
            #   w - vector of length n - optional weight of each phenotype
            #
            #   EFFECTS
            #   validates the inputs, then creates a GenomeScan object, storing the inputs in the 'Data' field
            #
            initialize = function(y, X, G, K, w = rep(1, length(y))) {

                #### CONDITIONS THAT CAUSE AN ERROR ####
                if (any(missing(y), missing(G))) {
                    stop('Must provide y (n-vector of phenotypes) and G (n-by-p matrix of genotypes) to define a GenomeScan.')
                }
                if (missing(K)) {
                    stop('No covariance information provided.  If none known, use lme4::lmer().  If known, input it as K.')
                }
                if (!all(length(y) == c(nrow(X), nrow(G), dim(K), length(w)))) {
                    stop("Input dimensions don't match.")
                }

                #### CONDITIONS THAT CAN BE MASSAGED INTO WORKING ####
                if (missing(X)) {
                    X <- matrix(data = 1, nrow = length(y))
                }

                # deal with NA data and non-positive weights
                to_carve <- is.na(y)
                to_carve <- to_carve | rowSums(x = is.na(X))
                to_carve <- to_carve | rowSums(x = is.na(G))
                to_carve <- to_carve | rowSums(x = is.na(K))
                to_carve <- to_carve | is.na(w)
                to_carve <- to_carve | w <= 0
                if (any(to_carve)) {
                    message('Removing ', sum(to_carve), ' observations due to NA phenotype, K, or weight or non-positive weight.')
                    y <- y[!to_carve]
                    X <- X[!to_carve,]
                    G <- G[!to_carve,]
                    K <- K[!to_carve, !to_carve]
                    w <- w[!to_carve]
                }

                Data <<- list(num_obs = length(y),
                              num_loci = ncol(G),
                              y = y, X = X, G = G, K = K, w = w)

                Results <<- list(LR = NULL, h2 = NULL)
            },


            #   INPUTS
            #   none
            #   TODO: allow method = eigen or cholesky
            #
            #   EFFECTS
            #   Does the computations that have to be done once per scan.
            #   Currently just eigen decomposes L
            #
            prep_scan = function() {

                message('Preparing GenomeScan...')
                L <- diag(1/sqrt(Data$w)) %*% Data$K %*% diag(1/sqrt(Data$w))
                eigen_L <- eigen(x = L, symmetric = TRUE)
                eigen_L$values <- check_eigen_decomposition(eigen_L)


                # fit the null model
                Intermediates_per_scan <<- list(L = L, eigen_L = eigen_L)
                null_fit <- fit_locus(locus_idx = NULL)
                Intermediates_per_scan <<- c(Intermediates_per_scan, null_LL = null_fit$LL)

                return(.self)
            },



            #   INPUTS
            #   none
            #   TODO: allow user to specify subset of chromosomes or loci
            #
            #   EFFECTS
            #   Computes 'Results' at each locus
            #
            conduct_scan = function() {

                if ('uninitializedField' %in% class(Intermediates_per_scan)) {
                    .self$prep_scan()
                }

                for (locus_idx in 1:Data$num_loci) {

                    fit <- fit_locus(locus_idx = locus_idx)

                    Results <<- list(LR = replace(x = Results$LR, list = locus_idx, values = 2*(fit$LL - Intermediates_per_scan$null_LL)),
                                     h2 = replace(x = Results$h2, list = locus_idx, values = fit$h2))
                }

                return(.self)

            },



            fit_locus = function(locus_idx) {

                if (is.null(locus_idx)) {
                    Intermediates_per_locus <<- list(XG = Data$X)
                } else {
                    Intermediates_per_locus <<- list(XG = cbind(Data$G[,locus_idx], Data$X))
                }

                opt <- optimize(f = .self$fit_locus_given_h2,
                                lower = 0,
                                upper = 1,
                                maximum = TRUE)

                return(list(h2 = opt[[1]],
                            LL = opt[[2]]))

            },




            fit_locus_given_h2 = function(h2) {

                calc_multiplier_eigen(h2 = h2)

                fit <- lm.fit(x = Intermediates_per_fit$M %*% Intermediates_per_locus$XG,
                              y = Intermediates_per_fit$M %*% Data$y)

                n <- Data$num_obs
                sigma2_mle <- sum(fit$residuals^2)/n
                loglik_uncorr <- (-0.5 * n * log(2*pi)) - (0.5 * n * log(sigma2_mle)) - (0.5 * n)
                loglik <- loglik_uncorr - 0.5 * Intermediates_per_fit$LDV

                return(loglik)
            },


            #   INPUTS
            #   h2 - heritability
            #
            #   EFFECTS
            #
            calc_multiplier_eigen = function(h2) {

                if (missing(h2)) {
                    stop('Must provide `h2` to `calc_multiplier_eigen`.')
                }

                eigen_L <- Intermediates_per_scan$eigen_L

                w <- Data$w
                M <- (1/sqrt(h2*eigen_L$values + (1 - h2)) * (1/sqrt(w) * t(eigen_L$vectors)))
                LDV <- sum(log(w)) + sum(log(h2*eigen_L$values + (1 - h2)))

                Intermediates_per_fit <<- list(M = M,
                                               LDV = LDV,
                                               h2 = h2)

            }


        )
    )