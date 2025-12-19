#' @name sam-adaptivelasso
#' 
#' @title Adaptive Lasso for Spatial Auto-regressive Model
#' @description Adaptive lasso for spatial auto-regressive model.
#' 
#' @param X Covariate matrix.
#' @param Y Response vector.
#' @param W Weight matrix.
#' @param rho Auto-correlation coefficient.
#' @param criterion *BIC* or *EBIC*.
#' @param lambda Tuning parameter.
#' @param lambda.vec Vector of \eqn{\lambda}.
#' @param max.iter Maximal number of iterations.
#' @param tol Covergence tolerance.
#' 
#' @return
#' * `sar_adaptivelasso`: `list(beta, sig2, rho, flag, BIC)`.
#' * `tune_sar_adaptivelasso`: `list(beta, sig2, rho, flag, BIC, lambda.opt)`.
#' 
#' @examples
#' set.seed(2025)
#' b0 <- c(1.5, 3.0, 2.0, rep(0.0, 3))
#' rho0 <- 0.2
#' sig0 <- 1.0
#' n <- 81
#' 
#' DF <- simu_sar_data_rook(b0, rho0, sig0, n)
#' Y <- DF[["y"]]
#' X <- as.matrix(DF[["X"]])
#' W0 <- DF[["W0"]]
#' 
#' system.time( rho.hat <- get_rho(X, Y, W0) )
#' system.time( tune_sar_adaptivelasso(X, Y, W0) )
#' system.time( tune_sar_adaptivelasso(X, Y, W0, rho.hat) )
#' system.time( tune_sar_adaptivelasso(X, Y, W0, criterion="EBIC") )
#' 
NULL
#> NULL


#' @rdname sam-adaptivelasso
#' @order 1
#' @export
sar_adaptivelasso <- function(X, Y, W, rho, lambda, max.iter = 1000, tol = 1e-6) {
    stopifnot( NROW(X) == length(Y) )
    stopifnot( length(lambda) == 1 & is.numeric(lambda) )
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    n <- NROW(X)
    p <- NCOL(X)
    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( length(Y) == n )

    A.rho <- diag(n) - rho * W
    y <- A.rho %*% Y
    
    ## initial estimation & active flag
    beta.ini <- tryCatch(
        coef(lm.fit(X, y, singular.ok = FALSE)),
        error = function(e) chol2inv(chol(crossprod(X) + diag(1e-3, p))) %*% crossprod(X, y)
    )

    beta_k <-beta.ini
    flag <- rep(TRUE, p)

    for (k in seq_len(max.iter)) {
        flag <- flag & (abs(beta_k) >= tol)

        lambda.tilde <- lambda / abs(beta.ini[flag])  # \tilde{λj} = λ / |βj0|, adaptive lasso
        D_k <- diag(c(lambda.tilde / abs(beta_k[flag])), length(lambda.tilde))
        
        beta.hat <- rep(0.0, p)
        beta.hat[flag] <- solve(
            crossprod(X[, flag, drop=FALSE]) + D_k,
            crossprod(X[, flag, drop=FALSE], y)
        )
        
        if (max(abs(beta.hat - beta_k)) < tol)
            break
        
        beta_k <- beta.hat
    }

    sig2.hat <- mean( (A.rho %*% Y - X[, flag, drop=FALSE] %*% beta.hat[flag])^2 )
    BIC <- n * log(sig2.hat) - 2*log(det(A.rho)) + log(n) * sum(flag)
    ebic.r <- max(0.0, log(p) / (2.0*log(n)))
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, sum(flag))

    egg <- list(
        beta = beta.hat,
        sig2 = sig2.hat,
        rho = rho,
        flag = flag,
        BIC = BIC,
        EBIC = EBIC
    )

    class(egg) <- "sam.adaptivelasso"
    return(egg)
}



#' @rdname sam-adaptivelasso
#' @order 2
#' @export
tune_sar_adaptivelasso <- function(
    X, Y, W, rho, criterion="BIC",
    lambda.vec = 10^seq(-1, 1, length = 100),
    max.iter = 1000,
    tol = 1e-6) {

    rho.hat <- ifelse(missing(rho), get_rho(X, Y, W), rho)

    idx <- sapply(lambda.vec, function(lambda)
        get(criterion, sar_adaptivelasso(X, Y, W, rho.hat, lambda, max.iter, tol))
    ) |> which.min()
    
    lambda.opt <- lambda.vec[idx]
    egg <- sar_adaptivelasso(X, Y, W, rho.hat, lambda.opt, max.iter, tol)
    egg["lambda.opt"] <- lambda.opt

    class(egg) <- "tune.sam.adaptivelasso"
    return(egg)
}
