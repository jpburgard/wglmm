#' @title Optimize likelihood of lme4 for transformation parameter
#' @name .optim_lmer
#' @keywords internal
#' @import lme4
#' @rawNamespace import(Matrix, except = c(cov2cor, toeplitz, update) )
#' @import stats
.optim_lmer <- function( formula, data = NULL, transform, ...){
  # Data Editing to be done by lme4
  mc <- match.call()
  mf <- match.call(expand.dots = TRUE)
  m <- match(c("data", "weights", "subset", "na.action", "contrasts", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::get_all_vars)
  mf$formula <- subbars(formula)
  dat <- eval(mf, parent.frame())
  rownames(dat) <- 1:nrow(dat)

  nam <- c("data", "weights", "subset", "na.action", "contrasts", "offset")[c(TRUE, TRUE, m[-(1:2)] != 0)]

  if( m[2] == 0 || length(unique(dat$weights)) == 1) w <- NULL else w <- dat$weights
  if(length(w) > 1) w <- w/sum(w)*length(w)

  nam <- c( "formula", nam)
  nam2 <- c( "formula", "dat", "w", nam[-(1:3)])
  env <- lapply(nam2, function(x) get(x))
  names(env) <- nam
  if(!is.null(w)) env$weights <- w
  env$formula <- formula
  # env$REML <- FALSE
  mod <- do.call(lmer, env)
  x <- 1
  form <- as.character(formula[2])
  if(transform != "none"){
    if(min(mod@frame[,form]) <= 0){
      mod@frame[,form] <- mod@frame[,form] - min(mod@frame[,form]) + 0.01
    }
    ys <- mod@frame[,form]
  }
  if(is.null(w)) w <- 1

  if(transform == "dual"){
    fu <- function(x){
      forms <- paste0("I( (", form, "^(", x, ") - 1/", form, "^(", x, ") )/(2 * (", x, ")) ) ~.")
      if(x == 0) forms <- paste0("log( ", form, " ) ~ .")
      mod <- update(mod, as.formula(forms))
      lls <- logLik(mod) + sum( log(ys^(x-1) + ys^(-x-1) )  * w )
      return(lls)
    }
    opt <- optimize(fu, lower = 0, upper = 10, maximum = TRUE)
    x <- opt$maximum
    forms <- paste0("I( (", form, "^(", x, ") - 1/", form, "^(", x, ") )/(2 * (", x, ")) ) ~.")
    if(x == 0) forms <- paste0("log( ", form, " ) ~ .")
    mod <- update(mod, as.formula(forms))
  
    }else if(transform == "box-cox"){
    fu <- function(x){
      forms <- paste0("I( (", form, "^(", x, ") - 1)/(", x, ")) ~.")
      if(x == 0) forms <- paste0("log( ", form, " ) ~ .")
      mod <- update(mod, as.formula(forms))
      lls <- logLik(mod) + (x-1) * sum( log(ys)  * w )
      return(lls)
    }
    opt <- optimize(fu, lower = -10, upper = 10, maximum = TRUE)
    x <- opt$maximum
    forms <- paste0("I( (", form, "^(", x, ") - 1)/(", x, ")) ~.")
    if(x == 0) forms <- paste0("log( ", form, " ) ~ .")
    mod <- update(mod, as.formula(forms))
  
    }
  return(list(transform = transform, lambda = x, model = mod))
}
