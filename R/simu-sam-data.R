#' @title Simulate Data
#' 
#' @description Rook weight matrix.
#' 
#' @param b0 True vector of \eqn{\beta_0}.
#' * `c(1.5, 3.0, 2.0, rep(0.0, 3))`
#' * `c(rep(1.5, 3), rep(2.0, 4), rep(3.0, 3), rep(0.0, 10))`
#' * `c(rep(1.5, 5), rep(2.0, 5), rep(3.0, 5), rep(0.0, 15))`
#' @param rho0 \eqn{\rho_0}.
#' @param sig0 \eqn{\sigma_0} for \eqn{\varepsilon}.
#' @param n Sample size.
#' 
#' @export
simu_sam_data_rook <- function(b0, rho0, sig0, n) {

    p <- length(b0)
    X <- data.frame( matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-"))) )

    W0 <- matrix(0.0, n, n)
    diag(W0[2:n, 1:(n-1)]) <- diag(W0[1:(n-1), 2:n]) <- 1.0

    y <- solve(diag(n) - rho0*W0, as.matrix(X) %*% b0 + rnorm(n, sd=sig0))

    return( list(y=y, X=X, W0=W0, s0=names(X)[which(b0 != 0)]) )
}
