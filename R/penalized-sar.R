#' @name penalized.sar
#' 
#' @title Penalized Methods for Spatial Auto-Regressive Model
#' 
#' @description
#' The penalized methods for the SAR model is developed for numerical simulations.
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix (row-sum scaled being one).
#' @param rho Auto-correlation coefficient.
#' @param lambda Tuning parameter.
#' @param lambda.vec Vector of \eqn{\lambda}.
#' @param max.iter Maximal number of iterations.
#' @param tol Covergence tolerance.
#' 
#' @return `list(beta, sig2, rho, lambda, BIC, EBIC)` for `bic / ebic`.
#' 
#' @examples
#' set.seed(2025)
#' b0 <- c(1.5, 3.0, 2.0, rep(0.0, 3))
#' rho0 <- 0.2
#' sig0 <- 1.0
#' n <- 81
#' 
#' DF <- simu_sar_data_rook(b0, rho0, sig0, n, sqrt(n), sqrt(n))
#' y <- DF[["y"]]
#' x <- DF[["x"]]
#' w <- DF[["w"]]
#' 
#' 
#' system.time( rho.hat <- get_rho(x, y, w) )
#' y.tilde <- (diag(NROW(x)) - rho.hat* w) %*% y
#' 
#' ## ---------- adaptive lasso ----------
#' system.time( tune_sar_adaptivelasso(x, y, w) )
#' system.time( egg <- tune_sar_adaptivelasso(x, y, w, rho.hat) )
#' 
#' flag <- egg[["bic.flag"]]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% egg[["bic.beta"]][flag]))^2 )
#' print( egg[["bic.sig2"]] - sig2.hat )
#' 
#' flag <- egg[["ebic.flag"]]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% egg[["ebic.beta"]][flag]))^2 )
#' print( egg[["ebic.sig2"]] - sig2.hat )
#' 
#' ## ---------- scad ----------
#' system.time( tune_sar_scad(x, y, w) )
#' system.time( tune_sar_scad(x, y, w, rho.hat) )
#' system.time( tune_sar_scad(x, y, w, lambda.vec=10^seq(-1, 1, length=5)) )
#' system.time( egg <- tune_sar_scad(x, y, w, rho.hat, lambda.vec=10^seq(-1, 1, length=5)) )
#' 
#' flag <- egg[["bic.flag"]]
#' y.hat <- drop(x[, flag, drop=FALSE] %*% egg[["bic.beta"]][flag]) + egg[["bic.intercept"]]
#' sig2.hat <- mean( (y.tilde - y.hat)^2 )
#' print( egg[["bic.sig2"]] - sig2.hat )
#' 
#' flag <- egg[["ebic.flag"]]
#' y.hat <- drop(x[, flag, drop=FALSE] %*% egg[["ebic.beta"]][flag]) + egg[["ebic.intercept"]]
#' sig2.hat <- mean( (y.tilde - y.hat)^2 )
#' print( egg[["ebic.sig2"]] - sig2.hat )
#' 
#' ## ---------- lasso ----------
#' system.time( tune_sar_lasso(x, y, w) )
#' system.time( tune_sar_lasso(x, y, w, rho.hat) )
#' system.time( tune_sar_lasso(x, y, w, lambda.vec=10^seq(-1, 1, length=5)) )
#' system.time( egg <- tune_sar_lasso(x, y, w, rho.hat, lambda.vec=10^seq(-1, 1, length=5)) )
#' 
#' flag <- egg[["bic.flag"]]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% egg[["bic.beta"]][flag]))^2 )
#' print( egg[["bic.sig2"]] - sig2.hat )
#' 
#' flag <- egg[["ebic.flag"]]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% egg[["ebic.beta"]][flag]))^2 )
#' print( egg[["ebic.sig2"]] - sig2.hat )
#' 
NULL
#> NULL
