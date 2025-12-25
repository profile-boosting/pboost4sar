#' @rdname penalized.sar
#' @export
tune_sar_lasso <- function(x, y, w, rho, lambda.vec,
                           standardize = FALSE, intercept = FALSE) {
    n <- length(y)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    rho <- ifelse(missing(rho), get_rho(x, y, w), rho)
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    A.rho <- diag(n) - rho * w
    y <- A.rho %*% y
    stopifnot( length(y) == n )

    if (missing(lambda.vec))
        models <- glmnet(x, y, "gaussian", standardize = standardize, intercept = intercept)
    else
        models <- glmnet(x, y, "gaussian", standardize = standardize, intercept = intercept, lambda = lambda.vec)

    lambda <- models[["lambda"]]
    beta.hat <- models[["beta"]]

    sig2.hat <- (as.numeric(y) - as.matrix(x %*% beta.hat))^2 |> colMeans() |> as.numeric()
    stopifnot( length(sig2.hat) == length(lambda) )

    dof <- models[["df"]]
    stopifnot( length(dof) == length(lambda) )

    BIC <- n * log(sig2.hat) + log(n) * dof # - 2*log(det(A.rho))
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, dof)

    opt.bic <- which.min(BIC)
    opt.ebic <- which.min(EBIC)
    egg <- list(
        bic.beta = as.numeric(beta.hat[, opt.bic]),
        bic.sig2 = sig2.hat[opt.bic],
        bic.rho = rho,
        bic.lambda = lambda[opt.bic],
        bic.BIC = BIC[opt.bic],
        bic.EBIC = EBIC[opt.bic],
        # --------------------
        ebic.beta = as.numeric(beta.hat[, opt.ebic]),
        ebic.sig2 = sig2.hat[opt.ebic],
        ebic.rho = rho,
        ebic.lambda = lambda[opt.ebic],
        ebic.BIC = BIC[opt.ebic],
        ebic.EBIC = EBIC[opt.ebic]
    )

    class(egg) <- "tune.sar.lasso"
    return(egg)
}