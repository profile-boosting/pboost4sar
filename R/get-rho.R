#' @title MLE of \eqn{\rho}
#' 
#' @description MLE of \eqn{\rho}.
#' 
#' @param x Covariate matrix.
#' @param y Response vector.
#' @param w Weight matrix (row-sum scaled).
#' 
#' @return MLE of \eqn{\rho}.
#' 
#' @examples
#' set.seed(2025)
#' b0 <- c(1.5, 3.0, 2.0, rep(0.0, 3))
#' rho0 <- 0.2
#' sig0 <- 1.0
#' n <- 81
#' 
#' DF <- simu_sar_data_rook(b0, rho0, sig0, n, sqrt(n), sqrt(n))
#' x <- DF[["x"]]
#' y <- DF[["y"]]
#' w <- DF[["w"]]
#' 
#' system.time( get_rho(x, y, w) )
#' 
#' @export
get_rho <- function(x, y, w) {
    stopifnot( is.matrix(x) )
    stopifnot( NROW(x) == length(y) )
    stopifnot( all( abs(rowSums(w) - 1.0) <= 1e-12 ) )
    stopifnot( all(w >= 0.0) )

    if (NCOL(x) > NROW(x))
        warning("Singularity Possible!")

    minusloglik <- function(rho) {
        A.rho <- diag(NROW(x)) - rho * w
        sig2.hat <- mean( residuals(lm.fit(x, A.rho %*% y))^2 ) # singular.ok=FALSE
        0.5*NROW(x) * log(sig2.hat) - log(det(A.rho))
    }

    optimize(minusloglik, c(-1, 1))$minimum
}