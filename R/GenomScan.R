setRefClass(Class = GenomeScan,
            fields = c('Data',            # y, X, K, and weights -- things that have interpretable meaning
                       'Intermediates',   # L, eigen_L, M, etc -- uninterpretable computational tools
                       'Results'),        # effect estimates, SEs, and p-values
            methods = list(

                initialize = function(y, X, K, weights = rep(1, length(y))) {

                    if (any(missing(y), missing(X))) { stop('Must provide y and X to do regression.') }
                    if (missing(K)) { stop('No covariance information provided.  Use lme4::lmer()') }

                    if (!all(length(y) == nrow(X), length(y) == dim(K), length(y) == length(weights))) {
                        stop("Input dimensions don't match.")
                    }


                },

                conduct_scan = function() {

                }

            )
)