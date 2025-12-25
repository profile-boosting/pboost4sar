#' @name simu-data
#' @title Simulate Data
#' 
#' @description
#' - `simu_sar_data_rook`: Rook weight matrix.
#' - `simu_sar_data_case`: Case weight matrix.
#' 
#' @param b0 True vector of \eqn{\beta_0}.
#' * `c(1.5, 3.0, 2.0, rep(0.0, 3))`
#' * `c(rep(1.5, 3), rep(2.0, 4), rep(3.0, 3), rep(0.0, 10))`
#' * `c(rep(1.5, 5), rep(2.0, 5), rep(3.0, 5), rep(0.0, 15))`
#' @param rho0 \eqn{\rho_0}.
#' @param sig0 \eqn{\sigma_0} for \eqn{\varepsilon}.
#' @param n Sample size.
#' 
#' @return `list(y, x, w)`.
#' 
NULL
#> NULL


#' @rdname simu-data
#' @order 1
#' @export
simu_sar_data_rook <- function(b0, rho0, sig0, n) {

    p <- length(b0)
    x <- matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-")))

    ## rNum: Number of row units
    ## cNum: Number of column units
    set_rook_matrix <- function(rNum, cNum=rNum){
        idxR <- as.vector( row(matrix(NA, rNum, cNum)) )
        idxC <- as.vector( col(matrix(NA, rNum, cNum)) )

        flag.row <- ( outer(idxR, idxR, \(i, j) i == j) & outer(idxC, idxC, \(i, j) abs(i - j) == 1) )
        flag.col <- ( outer(idxC, idxC, \(i, j) i == j) & outer(idxR, idxR, \(i, j) abs(i - j) == 1) )

        (flag.row | flag.col) * 1.0
    }
    W0 <- set_rook_matrix(sqrt(n))
    W0 <- W0 / rowSums(W0)
    stopifnot( NROW(W0) == n )
    stopifnot( NCOL(W0) == n )

    y <- solve(diag(n) - rho0*W0, x %*% b0 + rnorm(n, sd=sig0)) |> drop()

    return( list(y=y, x=x, w=W0) )
}



#' @rdname simu-data
#' @order 2
#' @export
simu_sar_data_case <- function(b0, rho0, sig0, n) {

    p <- length(b0)
    x <- matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-")))

    # W = I_R \otimes B_m
    R <- 10
    m <- 20
    W0 <- kronecker(diag(R), (matrix(1.0, m, m) - diag(m)) / (m - 1.0))
    W0 <- W0 / rowSums(W0)
    stopifnot( NROW(W0) == n )
    stopifnot( NCOL(W0) == n )

    y <- solve(diag(n) - rho0*W0, x %*% b0 + rnorm(n, sd=sig0)) |> drop()

    return( list(y=y, x=x, w=W0) )

}
