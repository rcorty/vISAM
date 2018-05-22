#' @title check_eigen_decomposition
#'
#' @param e the eigen decomposition to check
#' @param tol the threshold below which a number is said to be effectively zero, defaults to 1e-6
#'
#' @description Grabbed from MASS.  Useful to sparsify matrices when some eigenvalues are essentially zero.
#'
#' @return The eigen values with any values with absolute value less than tol zeroed.
#' @export
#'
check_eigen_decomposition <- function(e, tol = 1e-6) {

    # after MASS function mvrnorm
    if (!all(e$values >= -tol * abs(e$values[1L]))){
        stop("K is not positive definite")
    }
    if (any(e$values < 0)) {
        if (any(e$values < -tol)) {
            message("Zeroing negative eigenvalues: smallest eigenvalue was ",
                    min(e$values), "\n")
        }
        e$values <- pmax(e$values, 0)
    }
    return(e$values)
}


low_rank_eigen <- function(Z, K) {
    e <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE)
    return(list(values=e$values, vectors=qr.Q(qr(Z%*%e$vectors))))
}




subset_eigen <- function(l, good_idxs) {

    return(list(values = l$values[good_idxs],
                vectors = l$vectors[good_idxs, good_idxs]))
}