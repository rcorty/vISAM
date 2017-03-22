#' @title GenomeScan
#'
#' @description A Reference Class implementing a Genome Scan
#'
#' @field Data Things the user inputs.
#' They have interpretable meaning and define the GenomeScan.
#' Currently: y, X, G, K, and weights.
#' @field Intermediates_per_scan Things the GenomeScan will compute once per scan.
#' They are mathetmatical tools that can't really be interpreted.
#' Currently: L, eigen_L, and LL_null.
#' @field Intermediates_per_locus Things the GenomeScan will compute once per locus.
#' They are mathetmatical tools that can't really be interpreted.
#' Currently: XG
#' @field Intermediates_per_fit Things the Genomecan will compute many times per locus (once per trial fit on that locus).
#' These are interpretable but rapidly changing and not guaranteed to be finalized or optimal.
#' Currently: M, LDV, and h2
#' @field Results The results of the GenomeScan.
#' Currently: The h2 that maximizes the LL at each locus and the LR as compared with the no-locus (null) model.
#'
#'
#' @export
#'
GenomeScan <- setRefClass(
        Class = 'GenomeScan',
        fields = c('Data',                    # y, X, G, K, and weights --
                   'Intermediates_per_scan',  # L, eigen_L, LL_0 -- uninterpretable computational tools that are constant for a GenomeScan
                   'Intermediates_per_locus', # M -- uninterpretable computational tools that are constant for a locus
                   'Intermediates_per_fit',   # sigma_a and sigma_e -- interpretable, but changing very fast, don't look here for results
                   'Results'))                # maximum likelihood LL, sigma_a, and sigma_e for each locus



#' @title initialize
#'
#' @name GenomeScan_initialize
#'
#' @description Initialize a GenomeScan
#'
#' @param y vector of length n - the phenotype of each of n individuals
#' @param X matrix of dimension n-by-c - the covariate value of each individual for c covariates
#' @param G matrix of dimension n-by-p - the genotype of each individual at p loci
#' @param K matrix of dimension n-by-n - the covariance of the phenotype
#' @param w matrix of dimension n-by-n - the covariance of the phenotype
#'
#' @return an object of class GenomeScan
NULL
GenomeScan$methods(
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
    }
)


#' @title prep_scan
#'
#' @name GenomeScan_prep_scan
#'
#' @description Prepare a GenomeScan for running.  Does all the computations that need to be done exactly once per genome scan.
#'
#' @return an object of class GenomeScan
NULL
GenomeScan$methods(
    prep_scan = function() {

        message('Preparing GenomeScan...')
        L <- diag(1/sqrt(Data$w)) %*% Data$K %*% diag(1/sqrt(Data$w))
        eigen_L <- eigen(x = L, symmetric = TRUE)
        eigen_L$values <- check_eigen_decomposition(eigen_L)


        # fit the null model
        Intermediates_per_scan <<- list(L = L, eigen_L = eigen_L)
        null_fit <- .self$fit_locus(locus_idx = NULL)
        Intermediates_per_scan <<- c(Intermediates_per_scan, LL_null = null_fit$LL)

        return(.self)
    }
)



#' @title conduct_scan
#'
#' @name GenomeScan_conduct_scan
#'
#' @description Conducts the GenomeScan.
#'
#' @note TODO: allow user to specify subset of chromosomes or loci
#'
#' @return an object of class GenomeScan
#'
NULL
GenomeScan$methods(
    conduct_scan = function() {

        if ('uninitializedField' %in% class(Intermediates_per_scan)) {
            .self$prep_scan()
        }

        for (locus_idx in 1:Data$num_loci) {

            fit <- fit_locus(locus_idx = locus_idx)

            Results <<- list(LR = replace(x = Results$LR, list = locus_idx, values = 2*(fit$LL - Intermediates_per_scan$LL_null)),
                             h2 = replace(x = Results$h2, list = locus_idx, values = fit$h2))
        }

        return(.self)

    }
)




#' @title fit_locus
#'
#' @name GenomeScan_fit_locus
#'
#' @description Fit one locus of a GenomeScan.  Should not typically be called by a user.
#'
#' @return an object of class GenomeScan
#'
NULL
GenomeScan$methods(
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

    }
)


#' @title fit_locus_given_h2
#'
#' @name GenomeScan_fit_locus_given_h2
#'
#' @description Fit one locus at a specified value of h2.  Should not typically be called by a user.
#'
#' @return an object of class GenomeScan
#'
NULL
GenomeScan$methods(
    fit_locus_given_h2 = function(h2) {

        .self$calc_multiplier_eigen(h2 = h2)

        fit <- lm.fit(x = Intermediates_per_fit$M %*% Intermediates_per_locus$XG,
                      y = Intermediates_per_fit$M %*% Data$y)

        n <- Data$num_obs
        sigma2_mle <- sum(fit$residuals^2)/n
        loglik_uncorr <- (-0.5 * n * log(2*pi)) - (0.5 * n * log(sigma2_mle)) - (0.5 * n)
        loglik <- loglik_uncorr - 0.5 * Intermediates_per_fit$LDV

        return(loglik)
    }
)


#' @title calc_multiplier_eigen
#'
#' @name GenomeScan_calc_multiplier_eigen
#'
#' @description Compute a multiplier (aka rotation) matrix.  Details in in h2lmm_math_RWC.Rmd.
#'
#' @return an object of class GenomeScan
#'
NULL
GenomeScan$methods(
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
