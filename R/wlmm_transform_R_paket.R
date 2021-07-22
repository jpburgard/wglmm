#' @useDynLib wmm
#' @name wlmm_transform
#' @export
#' @title Survey-weighted Linear Mixed Models under Box-Cox and Dual Transformations
#' @aliases wlmm_transform
#' @description Estimates linear mixed models under Box-Cox and Dual Transformations with flexible random effects structure and survey weights.
#' @usage
#' wlmm_transform(formula, data = NULL, transform = "none",
#'        lambda = NULL, weights = NULL,
#'        iter1  = 501, iter2 = 501, MI = 500,
#'        tol1 = 3e-5, tol2 = 1e-5,
#'        trace = FALSE, nDecrease = 2, ...)
#' @param formula Formula for the Model to be estimated. Has the shape of lmer in the package lme4.
#' @param data Dataframe containing the variables
#' @param transform Chosen transformation: "none", "box-cox" or "dual"
#' @param lambda Start value for transformation parameter lambda
#' @param weights Survey weights
#' @param iter1 Maximum number of MCEM iterations
#' @param iter2 Maximum number of internal optimizations within a MCEM-step
#' @param MI Number of importance-sampled random numbers in E-step
#' @param tol1 Convergence criterion for MCEM algorithm
#' @param tol2 Convergence criterion for opimizations within a MCEM-step
#' @param trace Print output of each MCEM-step?
#' @param nDecrease How often is it allowed that the estimate of LL does not increase?
#' @param ... Additional options for model set-up
#' @details For comparability between runs, set a seed due to stochastic optimization.
#' @return A list with the following elements
#' \itemize{
#' \item model C++ output of estimation procedure:
#' \itemize{
#' \item coef Vector of fixed effects regression parameters
#' \item vcov_coef Fisher information of fixed effects parameters
#' \item lambda Transformation Parameter
#' \item VarCov List of random effects variance-covariance matrices
#' \item sigma Residual standard deviation
#' \item transform Chosen transformation
#' \item RE_mat Matrix with simulated random effects from last E-step, including importance weights
#' \item RE_mode Modes of the random effects
#' \item RE_invHess Inverse hessian of LL at mode of random effects
#' \item LLmod Joint maximum log-likelihood of observed data and the mode of random effects
#' \item LLexp Expected log-likelihood given the observed data
#' \item LL Estimated log-likelihood
#' \item convergence Has the model converged?
#' }
#' \item formula regression model
#' \item data Final estimate of transformation parameter
#' }
#' @references Burgard, Jans Pablo and Doerr, Patricia (2019). Data-driven transformations and survey-weighting for linear mixed models. Trier University. Working Paper.
#' @import lme4
#' @rawNamespace import(Matrix, except = c(cov2cor, toeplitz, update) )
#' @import stats
#' @examples
#' \dontrun{
#' library(emdi)
#' data("eusilcA_smp")
#' mod <- wlmm_transform( eqIncome ~ gender + eqsize + cash + self_empl +
#'                unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow +
#'                house_allow + cap_inv + tax_adj + (1|district), data = eusilcA_smp,
#'                transform = "box-cox", weights = eusilcA_smp$weight,
#'                MI = 500, iter1 = 250, trace = TRUE)
#' }
wlmm_transform <- function(formula, data = NULL, transform = "none", lambda = NULL, weights = NULL, iter1  = 501, iter2 = 501,
                     MI = 500, tol1 = 3e-5, tol2 = 1e-5, trace = FALSE, nDecrease = 2, ...){

    if(transform != "none"){
    # Data Editing to be done by lme4
    mc <- match.call()
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("data", "weights", "subset", "na.action", "contrasts", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m[1:2])]
    mf[[1L]] <- quote(stats::get_all_vars)
    mf$formula <- subbars(formula)
    dat <- eval(mf, parent.frame())
    rownames(dat) <- 1:nrow(dat)

    nam <- c("data", "weights", "subset", "na.action", "contrasts", "offset")[c(TRUE, TRUE, m[-(1:2)] != 0)]

    if( m[2] == 0 || length(unique(dat$weights)) == 1) w <- NULL else{
      w <- dat$weights
      w <- w/sum(w) * length(w)
    }
    nam <- c( "formula", nam)
    nam2 <- c( "formula", "dat", "w", nam[-(1:3)])
    env <- lapply(nam2[-3], function(x) get(x))
    names(env) <- nam[-3]

    mf$formula <- formula

    # weights falsch abgespeichert, aber wichtig fuer Reihenfolge
    env$formula <- formula
    # env$REML <- FALSE
    mod <- do.call(lmer, env)
    env$weights <- w
    conLin <- do.call(lFormula, env )
    form <- as.character(env$formula[2])

    if(transform != "none"){
      if(min(mod@frame[,form]) <= 0){
        mod@frame[,form] <- mod@frame[,form] - min(mod@frame[,form]) + 0.01
      }
    }

    Ztlist <- getME(mod, "Ztlist")
    flist <- getME(mod, "flist")

    nam <- sapply( strsplit( names(Ztlist), ".", fixed = TRUE), function(x) return(x[1]) )
    ni <- sapply(flist, function(x) length(unique(x)))
    ni <- ni[match(unique(nam), names(flist))]
    qvec <- mapply( function(z, ns) nrow(z)/ns, Ztlist, ni[match(nam, names(ni))])
    qvec <- tapply(qvec, nam, sum)
    qvec <- qvec[names(ni)]
    ni <- as.integer(ni)
    nam <- unique(nam)

    qvec <- as.integer(qvec)
    q <- sum(qvec)

    n.r <- sum(qvec*ni)

    ngrps <- length(qvec)

    if(is.null(w)) w <- as.numeric(1) else w <- conLin$fr$`(weights)`
    w <- w/sum(w)*length(w)

    Z <- getME(mod, "Z")
    X <- getME(mod, "X")
    ys <- getME(mod, "y")

    if(is.null(lambda)){

      env$weights <- NULL
      env$transform <- transform
      opt <- do.call(.optim_lmer, env)
      lambda <- opt$lambda
      mod <- opt$model
      }else{
        if(transform != "none"){
          if(lambda != 0){
            if(transform == "dual"){
              forms <- paste0("I( (", form, "^(", lambda, ") - ", form, "^( -", lambda, ") )/(2 * ", lambda, ") ) ~." )
            }else if(transform == "box-cox"){
              forms <- paste0("I( (", form, "^(", lambda, ") - 1 )/", lambda, ") ~.")
            }
          }else{
            forms <- paste0("log( ", form, " ) ~." )
          }
          mod <- update(mod, forms)
        }
      }
    beta <- fixef(mod)
    scale <- summary(mod)$sigma^2
    CovRE <- lapply(VarCorr(mod), as.matrix)

    modus <- unlist( lapply(ranef(mod), function(x) as.numeric(t(x))))
    out <- .wlmm_transform_cpp( beta, transform, X = X, y = ys, Z = Z, CovRE = CovRE, weights = w, qvec = qvec, ni = ni,
                               scale = scale, lambda = lambda, MI = MI, iter1 = iter1, iter2 = iter2, tol1 = tol1,
                               tol2 = tol2, trace = trace, nDecrease = nDecrease )

  return( list(model = out, formula = formula, data = dat) )
    }else{
      m <- match.call(expand.dots = TRUE)
      m <- m[names(m) !="transform"]
      m <- m[-1]
      
      mod <- do.call("wglmm", args = as.list(m))

      return(mod)
    }

}


