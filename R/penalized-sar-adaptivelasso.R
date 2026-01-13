#' @rdname penalized.sar
#' @export
sar_adaptivelasso <- function(x, y.tilde, lambda, beta.ini, max.iter = 1000, tol = 1e-6) {
    stopifnot( length(lambda) == 1 & is.numeric(lambda) )

    n <- NROW(x)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    y <- y.tilde
    beta.k <- beta.ini
    flag <- rep(TRUE, p)

    for (k in seq_len(max.iter)) {
        flag <- flag & (abs(beta.k) >= tol)

        lambda.tilde <- lambda / abs(beta.ini[flag])  # \tilde{λj} = λ / |βj0|, adaptive lasso
        D_k <- as.numeric(lambda.tilde / abs(beta.k[flag]))

        beta.hat <- rep(0.0, p)
        beta.hat[flag] <- solve_ridge(x[, flag, drop=FALSE], y, D_k)
        
        # D_k <- diag(c(lambda.tilde / abs(beta.k[flag])), length(lambda.tilde))
        # beta.hat[flag] <- solve(
        #     crossprod(x[, flag, drop = FALSE]) + D_k,
        #     crossprod(x[, flag, drop = FALSE], y)
        # )

        if (max(abs(beta.hat - beta.k)) < tol)
            break
        
        beta.k <- beta.hat
    }

    sig2.hat <- mean( (y - x[, flag, drop = FALSE] %*% beta.hat[flag])^2 )
    BIC <- n * log(sig2.hat) + log(n) * sum(flag) # - 2*log(det(A.rho))
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, sum(flag))

    egg <- list(
        beta = beta.hat,
        sig2 = sig2.hat,
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
    rho <- ifelse(missing(rho) || is.null(rho), get_rho(x, y, w), rho)

    y.tilde <- drop(y - rho * (w %*% y))
    stopifnot( length(y.tilde) == NROW(x) )

    if (NROW(x) > NCOL(x)) {
        beta.ini <- tryCatch(
            coef(lm.fit(x, y.tilde, singular.ok = FALSE)),
            error = function(e) solve_ridge(x, y.tilde, 1e-3)
        )
    } else {
        beta.ini <- solve_ridge(x, y.tilde, 1e-3)
    }
    stopifnot( length(beta.ini) == NCOL(x) )

    tmp <- sapply(lambda.vec, function(lambda)
        sar_adaptivelasso(x, y.tilde, lambda, beta.ini, max.iter, tol)[c("BIC", "EBIC", "lambda")]
    )

    lambda <- as.numeric(tmp["lambda", ])
    bic <- as.numeric(tmp["BIC", ])
    ebic <- as.numeric(tmp["EBIC", ])

    bic.opt.lambda <- lambda[which.min(bic)]
    ebic.opt.lambda <- lambda[which.min(ebic)]

    model.bic <- sar_adaptivelasso(x, y.tilde, bic.opt.lambda, beta.ini, max.iter, tol)
    model.ebic <- sar_adaptivelasso(x, y.tilde, ebic.opt.lambda, beta.ini, max.iter, tol)

    egg <- list(
        bic.beta = as.numeric(model.bic[["beta"]]),
        bic.sig2 = model.bic[["sig2"]],
        bic.rho = rho,
        bic.flag = model.bic[["flag"]],
        bic.lambda = model.bic[["lambda"]],
        bic.BIC = model.bic[["BIC"]],
        bic.EBIC = model.bic[["EBIC"]],
        # --------------------
        ebic.beta = as.numeric(model.ebic[["beta"]]),
        ebic.sig2 = model.ebic[["sig2"]],
        ebic.rho = rho,
        ebic.flag = model.ebic[["flag"]],
        ebic.lambda = model.ebic[["lambda"]],
        ebic.BIC = model.ebic[["BIC"]],
        ebic.EBIC = model.ebic[["EBIC"]]
    )

    class(egg) <- "tune.sar.adaptivelasso"
    return(egg)
}
