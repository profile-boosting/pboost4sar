#' @title Forward Stepwise Spatial Auto-Regressive Model
#' 
#' @description Forward stepwise strategy for spatial auto-regressive model.
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix (row-sum scaled being one).
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
#' system.time( egg <- forward_sar(x, y, w) )
#' y.tilde <- (diag(NROW(x)) - egg[["rho"]] * w) %*% y
#' 
#' flag <- egg[["flag"]]
#' beta.hat <- egg[["beta"]][sort(names(egg[["beta"]]))]
#' sig2.hat <- mean( (y.tilde - drop(x[, flag, drop=FALSE] %*% beta.hat))^2 )
#' print( egg[["sig2"]] - sig2.hat )
#' 
#' @export
forward_sar <- function(x, y, w) {
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

    # `fitsarlm` is taken from `psar`
    fitsarlm <- function(x, y, w) {
        stopifnot( NROW(x) == length(y) )
        stopifnot( is.matrix(x) )

        rho.hat  <- get_rho(x, y, w)
        y.tilde <- y - rho.hat * (w %*% y)
        beta.hat <- coef(lm.fit(x, y.tilde, singular.ok = FALSE))
        err      <- y.tilde - x %*% beta.hat
        sig2.hat <- mean( err^2 )
        BIC  <- 2 * attr(rho.hat, "minusloglik") + log(NROW(x)) * sum(NCOL(x))
        EBIC <- BIC + 2.0 * ebic.r * lchoose(p, NCOL(x))

        list(
            beta = beta.hat,
            sig2 = sig2.hat,
            rho = rho.hat,
            residuals = err,
            level = EBIC
        )
    }
    
    idx <- c()
    egg <- fitsarlm(x[, idx, drop=FALSE], y, w)
    ebic <- egg$level
    maxK <- floor(0.8*n)
    while (TRUE) {
        fCand <- setdiff(seq(p), idx)
        logLikCand <- rep(NA_real_, length(fCand))
        ebicCand <- rep(NA_real_, length(fCand))
        for (iC in seq_along(fCand)) {
            tmp.idx <- c(idx, fCand[iC])
            tmp.egg <- fitsarlm(x[, tmp.idx, drop=FALSE], y, w)
            logLikCand[iC] <- attr(tmp.egg$rho, "minusloglik")
            ebicCand[iC] <- tmp.egg$level
        }
        fBest <- fCand[which.min(logLikCand)]
        tmp.idx <- c(idx, fBest)
        stopifnot(setdiff(tmp.idx, idx) == fBest)
        if (length(tmp.idx) > maxK) break
        
        tmp.egg <- fitsarlm(x[, tmp.idx, drop=FALSE], y, w)
        tmp.ebic <- tmp.egg$level

        if (tmp.ebic >= ebic) break

        egg  <- tmp.egg
        ebic <- tmp.ebic
        idx  <- tmp.idx
    }

    flag <- rep(FALSE, p)
    flag[idx] <- TRUE
    egg[["flag"]] <- flag

    class(egg) <- "fs.sar"
    return(egg)
}