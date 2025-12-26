#' @name pvs4sar
#' 
#' @title Profiled Variable Selection for Spatial Auto-regressive Model
#' 
#' @description Profile boosting variable selection for spatial auto-regressive (SAR) model
#' \deqn{\bm{y} = \rho W \bm{y} + X \beta + \bm{\varepsilon},}
#' where \eqn{\bm{y}, \bm{\varepsilon} \in \mathbb{R}^n}, \eqn{X \in \mathbb{R}^{n \times p}},
#' \eqn{\rho \in (-1, 1)}.
#' 
#' `pvs4sar_lagsarlm()` equips the PVS procedure for [spatialreg::lagsarlm].
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix (row-sum scaled being one).
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
#' Also see [pboost::pboost].
#' @param keep Initial set of features that are included in model fitting.
#' **If `keep` is specified, it should also be fully included in the RHS of `formula`.**
#' Also see [pboost::pboost].
#' @param verbose Print the procedure path? Also see [pboost::pboost].
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
#' DF <- simu_sar_data_rook(b0, rho0, sig0, n, sqrt(n), sqrt(n))
#' x <- DF[["x"]]
#' y <- DF[["y"]]
#' w <- DF[["w"]]
#' 
#' system.time( egg <- pvs4sar(x, y, w, verbose=TRUE) )
#' y.tilde <- (diag(NROW(x)) - egg[["rho"]] * w) %*% y
#' 
#' flag <- egg[["flag"]]
#' beta.hat <- egg[["beta"]][sort(names(egg[["beta"]]))]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% beta.hat))^2 )
#' print( egg[["sig2"]] - sig2.hat )
#' 
#' system.time(
#'  pvs4sar_lagsarlm(y ~ ., data.frame(y, x), mat2listw(DF[["w"]], style="W"), verbose=TRUE)
#' )
#' 
NULL
#> NULL



#' @rdname pvs4sar
#' @export
pvs4sar <- function(x, y, w, maxK = NULL, keep = NULL, verbose = FALSE) {
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

    cnames <- colnames(x)
    if (is.null(cnames))
        cnames <- colnames(x) <- paste0("X", 1:p)

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
            beta = beta.hat,
            sig2 = sig2.hat,
            rho = rho.hat,
            residuals = err,
            level = EBIC
        )
    }

    showiter <- function(verbose, i.max, level=obj.level)
        if (verbose)
            if (missing(i.max))
                message(sprintf("Initial model with level=%.3f", level))
            else
                message(sprintf("Adding feature `%s': model level=%.3f", cnames[i.max], level))

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

    flag <- rep(FALSE, p)
    flag[idx] <- TRUE
    obj[["flag"]] <- flag

    class(obj) <- "pvs4sar"
    return(obj)
}




#' @rdname pvs4sar
#' @export
pvs4sar_lagsarlm <- function(
    formula, data = list(), listw, na.action, Durbin, type,
    method = "eigen", quiet = NULL, zero.policy = NULL, interval = NULL,
    tol.solve = .Machine$double.eps, trs = NULL, control = list(),
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( !missing(formula) )
    stopifnot( !missing(data) )

    cl <- match.call()

    sar_template <- cl
    sar_template$stopFun <- NULL
    sar_template$keep <- NULL
    sar_template$maxK <- NULL
    sar_template$verbose <- NULL
    sar_template[[1L]] <- quote(lagsarlm)

    required_paras <- c("data", "listw", "Durbin", "type")
    for (ipara in required_paras)
        if (!is.null(cl[[ipara]]))
            sar_template[[ipara]] <- eval(cl[[ipara]], envir = parent.frame())

    fitFun <- function(formula, data) {
        call <- sar_template
        call$formula <- formula
        call$data <- data
        return( eval(call, parent.frame()) )
    }

    pboost(formula, data, fitFun, residuals, stopFun,
           keep = keep, maxK = maxK, verbose = verbose)
}



#' @title Extended BIC for SAM
#' @description EBIC for spatial auto-regressive model.
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
