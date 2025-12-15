#' @name plagsarlm
#' 
#' @title Profile Boosting Variable Selection for Spatial Autoregressive Model
#' @description Profile boosting variable selection for spatial autoregressive model.
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix.
#' @param formula Parameters passed to [spatialreg::lagsarlm].
#' @param data Parameters passed to [spatialreg::lagsarlm].
#' @param listw Parameters passed to [spatialreg::lagsarlm].
#' @param na.action Parameters passed to [spatialreg::lagsarlm].
#' @param Durbin Parameters passed to [spatialreg::lagsarlm].
#' @param type Parameters passed to [spatialreg::lagsarlm].
#' @param method Parameters passed to [spatialreg::lagsarlm].
#' @param quiet Parameters passed to [spatialreg::lagsarlm].
#' @param zero.policy Parameters passed to [spatialreg::lagsarlm].
#' @param interval Parameters passed to [spatialreg::lagsarlm].
#' @param tol.solve Parameters passed to [spatialreg::lagsarlm].
#' @param trs Parameters passed to [spatialreg::lagsarlm].
#' @param control Parameters passed to [spatialreg::lagsarlm].
#' @param stopFun Stopping rule for profile boosting, which has the form
#'  stopFun(object) to evaluate the performance of model object returned
#'  by fitFun, such as [EBIC] or [BIC].
#' @param maxK Maximal number of identified features.
#' If `maxK` is specified, it will supress `stopFun`, saying that the
#' profile boosting continues until the procedure identifies `maxK` features.
#' The pre-specified features in `keep` are counted toward `maxK`.
#' @param keep Initial set of features that are included in model fitting.
#' **If `keep` is specified, it should also be fully included in the RHS of `formula`.**
#' @param verbose Print the procedure path?
#' 
#' @return Model object fitted on the selected features.
#' 
#' @examples
#' library(spdep)
#' library(spatialreg)
#' 
#' set.seed(2025)
#' b0 <- c(1.5, 3.0, 2.0, rep(0.0, 3))
#' rho0 <- 0.2
#' sig0 <- 1.0
#' n <- 81
#' 
#' DF <- simu_sam_data_rook(b0, rho0, sig0, n)
#' 
#' data <- with(DF, data.frame(y, X))
#' W0 <- mat2listw(DF[["W0"]], style="W")
#' system.time( egg1 <- plagsarlm(y ~ ., data, W0, verbose=TRUE) )
#' 
#' system.time( egg2 <- plagsarlm2(as.matrix(DF[["X"]]), DF[["y"]], DF[["W0"]], verbose=TRUE) )
#' 
#' 
NULL
#> NULL



#' @rdname plagsarlm
#' @order 1
#' @export
plagsarlm <- function(
    formula, data = list(), listw, na.action, Durbin, type,
    method = "eigen", quiet = NULL, zero.policy = NULL, interval = NULL,
    tol.solve = .Machine$double.eps, trs = NULL, control = list(),
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {

    # # --- solution A ---
    # lagsarlm_args <- names(formals(spatialreg::lagsarlm))
    # user_args <- as.list(environment())
    # args_for_lagsarlm <- user_args[intersect(names(user_args), lagsarlm_args)]


    # # --- solution B ---
    # cl <- match.call()
    # cl[[1L]] <- quote(lagsarlm)

    # call_env <- as.environment(as.list(environment()))

    # lagsarlmArgs <- names(formals(spatialreg::lagsarlm))
    # callArgs <- as.list(sam_template)[intersect(names(as.list(sam_template)), lagsarlmArgs)]

    # fitFun <- function(formula, data) {
    #     cl$formula <- formula
    #     cl$data <- data
    #     eval(cl, envir = call_env)
    # }

    # --- END ---

    cl <- match.call(expand.dots = TRUE)
    sam_template <- cl
    sam_template$stopFun <- NULL
    sam_template$keep <- NULL
    sam_template$maxK <- NULL
    sam_template$verbose <- NULL
    sam_template[[1L]] <- quote(lagsarlm)
    fitFun <- function(formula, data){
        call <- sam_template
        call$formula <- formula
        return( eval(call, parent.frame()) )
    }

    # scoreFun <- function(object) {
    #     stopifnot( inherits(object, "Sarlm") )

    #     n0 <- attr(logLik(object), "nobs")
    #     rho <- get(object, "rho")

    #     eta <- rep(coef(object)["(Intercept)"], n0)
    #     vnames <- names(coef(object))[-(1:2)]
    #     if (length(vnames) > 0)
    #         eta <- eta + as.matrix(data[, vnames, drop=FALSE]) %*% coef(object)[vnames]

    #     ## explicitly call `y` from `data`
    #     drop( (diag(n0) - rho * W0) %*% data$y ) - eta
    # }

    pboost(formula, data, fitFun, residuals, stopFun,
           keep=keep, maxK=maxK, verbose=verbose)
}




#' @title Extended BIC for SAM
#' @description EBIC for spatial autoregressive model.
#' 
#' @param object See [pboost::EBIC].
#' @param p See [pboost::EBIC].
#' @param p.keep See [pboost::EBIC].
#' @param ... See [pboost::EBIC].
#' 
#' @return EBIC value.
#' 
#' @export
EBIC.Sarlm <- function(object, p, p.keep, ...) {
    stopifnot( inherits(object, "Sarlm") )

    if (missing(p))
        p <- get("p", envir=parent.frame())
    if (missing(p.keep))
        p.keep <- get("p.keep", envir=parent.frame())

    dof <- attr(logLik(object), "df")
    n0 <- attr(logLik(object), "nobs")
    ebic.r <- max( 0.0, 1.0 - log(n0) / (2.0*log(p)) )
    ebic.penalty <- ifelse(
        ebic.r <= 0.0,
        0.0,
        2.0 * ebic.r * lchoose(p - p.keep, dof - p.keep)
    )

    stopifnot( is.finite(ebic.penalty) )
    BIC(object) + ebic.penalty
}





#' @rdname plagsarlm
#' @order 2
#' @export
plagsarlm2 <- function(x, y, w, maxK = NULL, keep = NULL, verbose = FALSE) {
    stopifnot( is.matrix(x) )
    stopifnot( all( abs(rowSums(w) - 1.0) <= 1e-12 ) )
    stopifnot( all(w >= 0.0) )
    stopifnot( !anyNA(x) )
    stopifnot( !anyNA(y) )
    stopifnot( !anyNA(w) )

    n <- NROW(x)
    p <- NCOL(x)
    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( length(y) == n )

    fitsarlm <- function(x, y, w) {
        stopifnot( NROW(x) == length(y) )
        stopifnot( is.matrix(x) )

        rho.hat  <- get_rho(x, y, w)
        A.rho    <- diag(NROW(x)) - rho.hat * w
        y.tilde  <- A.rho %*% y
        beta.hat <- coef(lm.fit(x, y.tilde, singular.ok = FALSE))
        err      <- y.tilde - x %*% beta.hat
        sig2.hat <- mean( err^2 )
        BIC      <- n * log(sig2.hat) - 2*log(det(A.rho)) + log(NROW(x)) * sum(NCOL(x))
        EBIC <- BIC + 2.0 * ebic.r * lchoose(p, NCOL(x))

        list(
            beta.hat = beta.hat,
            sig2.hat = sig2.hat,
            rho.hat = rho.hat,
            residuals = err,
            level = EBIC
        )
    }

    showiter <- function(verbose, i.max, level=obj.level)
        if (verbose)
            if (missing(i.max))
                message(sprintf("Initial model with level=%.3f", level))
            else
                message(sprintf("Adding feature `%i': model level=%.3f", i.max, level))

    if (!is.null(maxK)) {
        stopifnot(maxK < NCOL(x))
        maxK <- min( maxK+length(keep), n-1, p-1 )
    }

    idx <- keep
    obj <- fitsarlm(x[, idx, drop = FALSE], y, w)
    obj.level <- obj[["level"]]
    showiter(verbose, level=obj.level)

    while (TRUE) {
        i.max <- which.max(abs(crossprod(x, obj[["residuals"]])))
        idx.tmp <- unique(c(idx, i.max))
        obj.tmp <- fitsarlm(x[, idx.tmp, drop=FALSE], y, w)
        obj.tmp.level <- obj.tmp[["level"]]
        showiter(verbose, i.max, obj.tmp.level)

        if (!is.null(maxK)) {
            if (length(idx.tmp) > maxK) break
        } else {
            if (obj.tmp.level > obj.level) {
                break
            }
        }

        obj <- obj.tmp
        idx <- idx.tmp
        obj.level <- obj.tmp.level
    }

    return(obj)
}