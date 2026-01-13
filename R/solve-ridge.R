#' @title Ridge-Type Estimation
#' 
#' @description Solve the linear system \eqn{(X^\top X + D) \beta = y}.
#' Using Woodbury matrix identity
#' \eqn{(D + X^\top X)^{-1} = D^{-1} - D^{-1} X^\top (I_n + X D^{-1} X'^\top)^{-1} X^\top D^{-1}},
#' we have
#' \deqn{(D + X^\top X)^{-1} X^\top Y = D^{-1} X^\top (I_n + X D^{-1} X^\top)^{-1} y.}
#' Specifically,
#' \deqn{(X^\top X + \lambda I_p)^{-1} X^\top y = X^\top (XX^\top + \lambda I_n)^{-1} y.}
#' 
#' @param X Matrix of size \eqn{n \times p}.
#' @param y Vector of length \eqn{n}.
#' @param d Positive scalar or vector.
#' 
#' @return Vector of length \eqn{p}.
#' 
#' @examples
#' set.seed(2026)
#' n <- 100
#' p <- 500
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' 
#' system.time( b1 <- solve_ridge(X, y, 0.1) )
#' max(abs(b1 - solve(crossprod(X) + 0.1 * diag(p), crossprod(X, y))))
#' 
#' d <- runif(p, 0.01, 0.5)
#' system.time( b2 <- solve_ridge(X, y, d) )
#' max(abs(b2 - solve(crossprod(X) + diag(d), crossprod(X, y))))
#' 
#' @export
solve_ridge <- function(X, y, d) {
    if (!is.matrix(X))X <- as.matrix(X)
    n <- NROW(X)
    p <- NCOL(X)
    stopifnot( is.numeric(X) & is.numeric(y) )
    stopifnot( min(n, p) > 0 )
    stopifnot( length(y) == n )
    stopifnot( all(d > 0) )
    stopifnot( length(d) == 1 || length(d) == NCOL(X) )

    if (n > p) {
        A <- crossprod(X)
        diag(A) <- diag(A) + d
        cA <- chol(A)
        return(drop( backsolve(cA, backsolve(cA, crossprod(X, y), transpose=TRUE)) ))
    } else {
        if (length(d) == 1) {
            A <- tcrossprod(X)
            diag(A) <- diag(A) + d
            cA <- chol(A)
            return(drop( crossprod(X, backsolve(cA, backsolve(cA, y, transpose=TRUE))) ))
        } else if (length(d) == NCOL(X)) {
            d_inv_X <- sweep(X, 2, d, "/")
            A <- tcrossprod(X, d_inv_X)
            diag(A) <- diag(A) + 1.0
            cA <- chol(A)
            return(drop( crossprod(d_inv_X, backsolve(cA, backsolve(cA, y, transpose=TRUE))) ))
        } else {
            stop("length(d) is not wrong.")
        }
    }
}
