#' @rdname penalized.sar
#' @export
tune_sar_scad <- function(x, y, w, rho, lambda.vec) {
    n <- length(y)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    rho <- ifelse(missing(rho) || is.null(rho), get_rho(x, y, w), rho)
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    A.rho <- diag(n) - rho * w
    y <- drop(A.rho %*% y)
    stopifnot( length(y) == n )

    if (missing(lambda.vec))
        models <- ncvreg(x, y, "gaussian", "SCAD")
    else
        models <- ncvreg(x, y, "gaussian", "SCAD", lambda = lambda.vec)

    lambda <- models[["lambda"]]

    idx.beta <- grep("(Intercept)", rownames(models$beta), fixed = TRUE, invert = TRUE)
    stopifnot( length(idx.beta) == p )
    beta.hat <- models[["beta"]][idx.beta, ]

    sig2.hat <- models[["loss"]] / models[["n"]] # `loss` is SSE for `gaussian`
    stopifnot( length(sig2.hat) == length(lambda) )

    dof <- colSums( abs(beta.hat) >= 1e-4 ) |> as.numeric()
    stopifnot( length(dof) == length(lambda) )

    BIC <- n * log(sig2.hat) + log(n) * dof # - 2*log(det(A.rho))
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, dof)

    opt.bic <- which.min(BIC)
    opt.ebic <- which.min(EBIC)

    bic.flag <- abs(as.numeric(beta.hat[, opt.bic])) > 1e-6
    ebic.flag <- abs(as.numeric(beta.hat[, opt.ebic])) > 1e-6
    stopifnot( length(bic.flag) == p )
    stopifnot( length(ebic.flag) == p )

    egg <- list(
        bic.beta = as.numeric(beta.hat[, opt.bic]),
        bic.intercept = models[["beta"]]["(Intercept)", opt.bic],
        bic.sig2 = sig2.hat[opt.bic],
        bic.rho = rho,
        bic.flag = bic.flag,
        bic.lambda = lambda[opt.bic],
        bic.BIC = BIC[opt.bic],
        bic.EBIC = EBIC[opt.bic],
        # --------------------
        ebic.beta = as.numeric(beta.hat[, opt.ebic]),
        ebic.intercept = models[["beta"]]["(Intercept)", opt.ebic],
        ebic.sig2 = sig2.hat[opt.ebic],
        ebic.rho = rho,
        ebic.flag = ebic.flag,
        ebic.lambda = lambda[opt.ebic],
        ebic.BIC = BIC[opt.ebic],
        ebic.EBIC = EBIC[opt.ebic]
    )

    class(egg) <- "tune.sar.scad"
    return(egg)
}
