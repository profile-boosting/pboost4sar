#' @rdname penalized.sar
#' @export
sar_adaptivelasso <- function(x, y, w, rho, lambda, max.iter = 1000, tol = 1e-6) {
    stopifnot( length(lambda) == 1 & is.numeric(lambda) )
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    n <- NROW(x)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    A.rho <- diag(n) - rho * w
    y <- drop(A.rho %*% y)
    stopifnot( length(y) == n )
    
    ## initial estimation & active flag
    beta.ini <- tryCatch(
        coef(lm.fit(x, y, singular.ok = FALSE)),
        error = function(e) drop(chol2inv(chol(crossprod(x) + diag(1e-3, p))) %*% crossprod(x, y))
    )

    beta_k <- beta.ini
    flag <- rep(TRUE, p)

    for (k in seq_len(max.iter)) {
        flag <- flag & (abs(beta_k) >= tol)

        lambda.tilde <- lambda / abs(beta.ini[flag])  # \tilde{λj} = λ / |βj0|, adaptive lasso
        D_k <- diag(c(lambda.tilde / abs(beta_k[flag])), length(lambda.tilde))
        
        beta.hat <- rep(0.0, p)
        beta.hat[flag] <- solve(
            crossprod(x[, flag, drop = FALSE]) + D_k,
            crossprod(x[, flag, drop = FALSE], y)
        )
        
        if (max(abs(beta.hat - beta_k)) < tol)
            break
        
        beta_k <- beta.hat
    }

    sig2.hat <- mean( (y - x[, flag, drop = FALSE] %*% beta.hat[flag])^2 )
    BIC <- n * log(sig2.hat) - 2*log(det(A.rho)) + log(n) * sum(flag)
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, sum(flag))

    egg <- list(
        beta = beta.hat,
        sig2 = sig2.hat,
        rho = rho,
        flag = flag,
        BIC = BIC,
        EBIC = EBIC,
        lambda = lambda
    )

    class(egg) <- "sar.adaptivelasso"
    return(egg)
}



#' @rdname penalized.sar
#' @export
tune_sar_adaptivelasso <- function(x, y, w, rho, lambda.vec = 10^seq(-1, 1, length = 100)) {

    max.iter <- 1000
    tol <- 1e-6
    rho <- ifelse(missing(rho), get_rho(x, y, w), rho)

    tmp <- sapply(lambda.vec, function(lambda)
        sar_adaptivelasso(x, y, w, rho, lambda, max.iter, tol)[c("BIC", "EBIC", "lambda")]
    )

    lambda <- as.numeric(tmp["lambda", ])
    bic <- as.numeric(tmp["BIC", ])
    ebic <- as.numeric(tmp["EBIC", ])

    bic.opt.lambda <- lambda[which.min(bic)]
    ebic.opt.lambda <- lambda[which.min(ebic)]

    model.bic <- sar_adaptivelasso(x, y, w, rho, bic.opt.lambda, max.iter, tol)
    model.ebic <- sar_adaptivelasso(x, y, w, rho, ebic.opt.lambda, max.iter, tol)

    egg <- list(
        bic.beta = as.numeric(model.bic[["beta"]]),
        bic.sig2 = model.bic[["sig2"]],
        bic.rho = model.bic[["rho"]],
        bic.flag = model.bic[["flag"]],
        bic.lambda = model.bic[["lambda"]],
        bic.BIC = model.bic[["BIC"]],
        bic.EBIC = model.bic[["EBIC"]],
        # --------------------
        ebic.beta = as.numeric(model.ebic[["beta"]]),
        ebic.sig2 = model.ebic[["sig2"]],
        ebic.rho = model.ebic[["rho"]],
        ebic.flag = model.ebic[["flag"]],
        ebic.lambda = model.ebic[["lambda"]],
        ebic.BIC = model.ebic[["BIC"]],
        ebic.EBIC = model.ebic[["EBIC"]]
    )

    class(egg) <- "tune.sar.adaptivelasso"
    return(egg)
}
