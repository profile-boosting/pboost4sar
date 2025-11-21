#' @name simu-data
#' @title Simulate Data
#' 
#' @description 
#' - `simu_sam_data_rook`: Rook weight matrix.
#' - `simu_sam_data_case`: Case weight matrix.
#' 
#' @param b0 True vector of \eqn{\beta_0}.
#' * `c(1.5, 3.0, 2.0, rep(0.0, 3))`
#' * `c(rep(1.5, 3), rep(2.0, 4), rep(3.0, 3), rep(0.0, 10))`
#' * `c(rep(1.5, 5), rep(2.0, 5), rep(3.0, 5), rep(0.0, 15))`
#' @param rho0 \eqn{\rho_0}.
#' @param sig0 \eqn{\sigma_0} for \eqn{\varepsilon}.
#' @param n Sample size.
#' @param m In Case spatial weight matrix, \eqn{W = I_R \otimes B_m}.
#' @param R In Case spatial weight matrix, \eqn{W = I_R \otimes B_m}.
#' 
#' @return `list(y, X, W0, s0)`.
#' 
NULL
#> NULL


#' @rdname simu-data
#' @order 1
#' @export
simu_sam_data_rook <- function(b0, rho0, sig0, n) {

    p <- length(b0)
    X <- data.frame( matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-"))) )

    # W0 <- matrix(0.0, n, n)
    # diag(W0[2:n, 1:(n-1)]) <- diag(W0[1:(n-1), 2:n]) <- 1.0
    swmrook<-function(num){
        Wr<-matrix(0,ncol=num,nrow=num)
        up<-sqrt(num)
        dw<-num-sqrt(num)+1
        for(i in 1:num){
            if(i>up) Wr[i,i-up]=1##上边相邻的单元
            if(i<dw) Wr[i,i+up]=1##下面相邻的单元
            if(i%%up!=1) Wr[i,i-1]=1##左边相邻单元
            if(i%%up!=0) Wr[i,i+1]=1##右边相邻单元
        }
        return(Wr)
    }
    W0 <- swmrook(n)
    W0 <- W0 / rowSums(W0)

    y <- solve(diag(n) - rho0*W0, as.matrix(X) %*% b0 + rnorm(n, sd=sig0))

    return( list(y=y, X=X, W0=W0, s0=names(X)[which(b0 != 0)]) )
}



#' @rdname simu-data
#' @order 2
#' @export 
simu_sam_data_case <- function(b0, rho0, sig0, m, R=8) {

    n <- m * R
    p <- length(b0)
    X <- data.frame( matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-"))) )

    W0 <- kronecker(diag(R), (matrix(1.0, m, m) - diag(m)) / (m - 1.0))
    W0 <- W0 / rowSums(W0)

    y <- solve(diag(n) - rho0*W0, as.matrix(X) %*% b0 + rnorm(n, sd=sig0))

    return( list(y=y, X=X, W0=W0, s0=names(X)[which(b0 != 0)]) )

}
