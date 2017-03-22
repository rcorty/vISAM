check_eigen_decomposition <- function(e, tol = 1e-6) {

    # after MASS function mvrnorm
    if (!all(e$values >= -tol * abs(e$values[1L]))){
        stop("K is not positive definite")
    }
    if (any(e$values < 0)) {
        if (any(e$values < -tol)) {
            warning("Zeroing negative eigenvalues: smallest eigenvalue was ",
                    min(e$values), "\n")
        }
        e$values <- pmax(e$values, 0)
    }
    return(e$values)
}
