#' @useDynLib wmm
#' @title Survey-weighted Generalized Linear Mixed Models
#' @aliases wglmm
#' @description Estimates generalized linear mixed models with flexible random effects structure and survey weights.
#' @usage wglmm(formula, data = NULL, family = gaussian(), weights = NULL,
#'   iter1  = 501, iter2 = 501,  MI = 500, tol1 = 2e-4,
#'   tol2 = 1e-5, trace = TRUE, nDecrease = 2, ...)
#' @param formula Formula for the Model to be estimated. Has the shape of lmer in the package lme4.
#' @param data Dataframe containing the variables
#' @param family Distribution from the exponential family. You can choose between gaussian(), binomial(),  poisson(),  Gamma() and  exponential(). However, only  gaussian() and binomial() are already extensively tested.
#' @param weights Survey weights
#' @param iter1 Maximum number of MCEM iterations
#' @param iter2 Maximum number of internal optimizations within a MCEM-step
#' @param MI Number of importance-sampled random numbers in E-step
#' @param tol1 Convergence criterion for MCEM algorithm
#' @param tol2 Convergence criterion for optimizations within a MCEM-step
#' @param trace Print output of each MCEM-step?
#' @param nDecrease How often is it allowed that the estimate of LL does not increase?
#' @param ... Additional options for model set-up
#' @details For comparability between runs, set a seed due to stochastic optimization.
#' @export
#' @return A list with the following elements
#' \itemize{
#' \item coef Vector of fixed effects regression parameters
#' \item VarCov List of random effects variance-covariance matrices
#' \item scale Scale parameter for exponential family
#' \item RE_mat Matrix with simulated random effects from last E-step, including importance weights
#' \item RE_mode List with modes of the random effects
#' \item residuals List with different types of residuals
#' \item LLmod Joint maximum log-likelihood of observed data and the mode of random effects
#' \item LLexp Expected log-likelihood given the observed data
#' \item niter Number of iterations
#' \item convergence Has MCEM converged?
#' }
#' @references Burgard, Jan Pablo and Doerr, Patricia (2018). Survey-weighted Generalized Linear Mixed Models. Trier University. Working Paper.
#' @examples
#' \dontrun{
#' library(mvtnorm)
#' n            <- 200
#' beta         <- c(4, -2, -1)
#' covmat       <- matrix(0.5, ncol = 2, nrow = 2)
#' diag(covmat) <- c(0.7, 1.3)
#' g1           <- 10
#' g2           <- 20
#' X1           <- rnorm(n, mean = 2)
#' X2           <- rexp(n, rate = 1)
#'
#' group1 <- rep(1:g1, length.out = n)
#' group2 <- rep(1:g2, length.out = n)
#' re1    <- rnorm(g1, sd = 2)
#' re2    <- rmvnorm(g2, sigma = covmat)
#'
#' modX  <- model.matrix( ~ X1 + X2)
#' modZ1 <- model.matrix( ~ -1 + as.factor(group1))
#' modZ2A <- model.matrix( ~ -1 + as.factor(group2))
#' modZ2B <- model.matrix( ~ -1 + X1:as.factor(group2))
#'
#' eta    <- modX %*% beta + modZ1 %*% re1 +
#'           modZ2A %*% as.vector(re2[,1]) +
#'           modZ2B %*% as.vector(re2[,2])
#'
#' lin    <- eta + rnorm(n, sd = 2.3)
#' dfsamp <- data.frame(group1 = factor(group1, levels = 1:10),
#'                      group2 = factor(group2, levels = 1:20),
#'                      X1 = X1, X2 = X2)
#' dfsamp$lin <- lin
#' modLin <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#'             data = dfsamp, trace = TRUE, iter1 = 250,
#'             iter2 = 1001, MI = 2000)
#' }
#' @import lme4
#' @rawNamespace import(Matrix, except = c(cov2cor, toeplitz, update) )
#' @import stats
wglmm <- function(formula, data = NULL, family = gaussian(),
                  weights = NULL, iter1  = 501, iter2 = 501,  MI = 500,
                     tol1 = 2e-4, tol2 = 1e-5,
                     trace = TRUE, nDecrease = 2, ...){

  # Data Editing to be done by lme4
  mc <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data", "weights", "subset", "na.action", "offset", "contrasts"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::get_all_vars)
  mf$formula <- subbars(formula)
  dat <- eval(mf, parent.frame())
  rownames(dat) <- 1:nrow(dat)

  nam <- c("data", "weights", "subset", "na.action", "contrasts", "offset")[c(TRUE, TRUE, m[-(1:2)] != 0)]

  if( m[2] == 0 || length(unique(dat$weights)) == 1) w <- NULL else w <- dat$weights

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()

  nam <- c( "formula", nam)
  nam2 <- c( "formula", "dat", "w", nam[-(1:3)])
  env <- lapply(nam2, function(x) get(x))
  names(env) <- nam

  mf$formula <- formula

  if (isTRUE(all.equal(family, gaussian()))) {
    # weights stored as heteroskedasticity correction, but correct ordering of elements
    conLin <- do.call(lFormula, env )
    env$weights <- w
    env$formula <- nobars( formula )
    mod <- do.call(lm, env )
    scale <- summary(mod)$sigma^2
  } else {

    env$family <- family
    conLin <- do.call(glFormula, env )
    env$weights <- w
    env$formula <- nobars( formula )
    mod <- do.call(glm, env )
    scale <- scaleNeu <- 1
    if(family$family == "Gamma") scale <- 1/summary(mod)$dispersion
  }
  # beta <- fixef(mod)
  # CovREalt <- CovRE <- lapply(VarCorr(mod), as.matrix)
  beta <- coef(mod)
  rm(mod)

  nam <- sapply( strsplit( names(conLin$reTrms$Ztlist), "| ", fixed = TRUE), function(x) return(x[2]) )
  ni <- sapply(conLin$reTrms$flist, function(x) length(unique(x)))

  ni <- ni[match(nam, names(conLin$reTrms$flist))]
  ni <- as.integer(ni)

  qvec <- mapply( function(z, ns) nrow(z)/ns, conLin$reTrms$Ztlist, ni)
  qvec <- as.integer(qvec)
  q <- sum(qvec)
  if(family$family == "Gamma")  CovREalt <- CovRE <- lapply(qvec, function(q) diag(q)*0.01) else CovREalt <- CovRE <- lapply(qvec, function(q) diag(q)*0.1) 

  n.r <- sum(qvec*ni)

  ngrps <- length(qvec)

  if(is.null(w)) w <- as.numeric(1) else w <- as.numeric( conLin$fr$`(weights)` )
  w <- w/sum(w)*length(w)

  X <- as.matrix(conLin$X)
  y <- conLin$fr[,as.character(formula[2])]
  if(is.matrix(y) && ncol(y) > 1) y <- y[,1]/rowSums(y)
  if(is.matrix(y)) y <- y[,1]
  if(is.factor(y)) y <- as.numeric(y)-1
  Z <- t( conLin$reTrms$Zt )

  n <- length(y)
  p <- length(beta)

  i <- 0
  eps <- 1.1 * tol1
  MI2 <- MI0 <- MI
  mi <- round(MI*0.1)

  # Initialization:
  phi.alt <- phi.altalt <- beta

  modus <- rep(0,n.r)

  #Chol <- diag(n.r)/sqrt( sum((.GradientAll(modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale))^2) )
   Chol <- 0.1 * diag(n.r)
  out <- .GenModusNeu( modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = Chol, scale)
  Chol <- out$invChol
  LLalt <- .loglikTot(out$modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)
  LLallTimeBest <- LL <- LLalt
  phiallTimeBest <- phi.alt
  CovREallTimeBest <- CovRE
  scaleallTimeBest <- scale
  counter <- 0

  while(i < iter1 && max(eps) > tol1){
    k <- 0
    increaseLL <- FALSE
    while(k < 10 & !increaseLL){
      Us <- .ImportanceSampling(MI = MI2 + k*mi, modus = out$modus, Chol = 1.05*Chol, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)

      # Step 1 of optimization: Fixed effects
      outFix <- try( .NewCoefImp( Us[,-(n.r+1)], phi.alt, family$family, X, y, Z, weights = w, weightsImp = Us[,n.r+1], tol = tol2, MaxIt = iter2, scale) )
      if( class(outFix) == "try-error" || !outFix$convergence || any(is.infinite(outFix$phi)) || any(is.na(outFix$phi)) ) break else phi.neu <- outFix$phi

      # Step 2 of optimization: RE-covariance matrix
      outRE <- .NewCovREImp( Us[,-(n.r+1)], qvec, ni, weightsImp = Us[,n.r+1])

        if(family$family == "gaussian"){
          #out <- try( .GenModusNeu( out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = out$invChol, scale = scale) )
          #if(class(out)=="try-error") break
          
          Xphi <- as.numeric(X %*% phi.neu)
          #linpred <- as.numeric(X %*% phi.neu + Z %*% out$modus)
          #res <- y - linpred
          scaleNeu <- sum( apply( Us[,1:n.r], 1, function(x) sum(w * (y - Xphi - as.numeric(Z %*% x))^2)/n ) * Us[,n.r+1] )
          #if(length(w) > 1) scaleNeu <- sum(w*res^2)/sum(w) else scaleNeu <- mean(res^2)
        }
        if(family$family == "Gamma"){
          # linpred <- as.numeric(X %*% phi.neu + Z %*% out2$modus)
          # res <- y - 1/linpred
          # if(length(w) > 1) scaleNeu <- sum(w)/sum(w*res^2*(linpred)^2) else scaleNeu <- 1/mean(res^2*(linpred)^2)
          Xphi <- as.numeric(X %*% phi.neu)
          scaleNeu <- sum( apply( Us[,1:n.r], 1, function(x){
            mu <- 1/(Xphi + as.numeric(Z %*% x))
            mu[mu <= 0] <- min(mu[mu > 0])*0.5
            # out <- sum(w*log(mu/y))*2/n # Deviance method
            out <- sum( w * (y - mu)^2/mu^2)/n # Pearson method
            # out <- sum( w * (y - mu)^2)/n
            return(out)} ) * Us[,n.r+1])
          # scaleNeu <- 1/scaleNeu

        }

      Us <- .ImportanceSampling(MI = MI2 + k*mi, modus = out$modus, Chol = 1.05*Chol, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)

      LL <- .loglikTotImp(Us[,1:n.r], phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, Us[,n.r+1], scaleNeu)
      LLalt <- .loglikTotImp(Us[,1:n.r], phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, Us[,n.r+1], scale)
      increaseLL <- LL < LLalt

      k <- k+1
      if(i == 0) break;
    }

    # Chol <- diag(n.r)/sqrt( sum((.GradientAll(out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, scaleNeu))^2) )
    out2 <- try( .GenModusNeu( out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = Chol, scale = scale) )
    if(class(out2)=="try-error") break
    Chol <- out2$invChol

    MI2 <- MI2 + mi
    if(any(sapply(outRE, det)== 0)) break;

    if( increaseLL | i == 0){
      LLallTimeBest <- LL
      phiallTimeBest <- phi.neu
      CovREallTimeBest <- outRE
      scaleallTimeBest <- scaleNeu
    }else{
      counter <- counter + 1
    }
    if(counter == nDecrease){
      warning("Likelihood does not increase anymore!")
      break
    }

    out <- out2
    rm(out2)

    eps1 <- abs(phi.neu-phi.alt)/(abs(phi.alt) + 1)

    if( i > 0){
      eps2 <- unlist( mapply( function(mat1, mat2) as.vector(abs(mat1[lower.tri(mat1, diag = TRUE)]-mat2[lower.tri(mat1, diag = TRUE)])/(abs(mat2[lower.tri(mat2, diag = TRUE)]) + 1)), outRE, CovRE) )
      eps22 <- unlist( mapply( function(mat1, mat2) as.vector(abs(mat1[lower.tri(mat1, diag = TRUE)]-mat2[lower.tri(mat1, diag = TRUE)])/(abs(mat2[lower.tri(mat2, diag = TRUE)]) + 1)), outRE, CovREalt) )
      } else eps2 <- eps22 <- max(0.05, 1.1*tol1)
    eps <- c(eps1, eps2)

    if( MI2 %% 2 == 1 ) MI2 <- MI2 + 1

    i <- i+1

    if(trace){
      print(paste0("M-Step ", i, " of ", iter1, "  Maximum variable change in iteration: ", max(eps)*100, " %") )
      print(cbind(phi.alt, phi.neu) )
      print(mapply( function(mat1, mat2) cbind(mat1, mat2), CovRE, outRE, SIMPLIFY = FALSE))
      print(paste0("Increased LL: ", increaseLL))
    }

    LLalt <- LL
    scale <- scaleNeu
    phi.altalt <- phi.alt
    CovREalt <- CovRE
    CovRE <- outRE
    phi.alt <- phi.neu

  }

  # Fixed Effects Editing
  coefs <- as.numeric( phiallTimeBest )
  CovRE <- CovREallTimeBest
  scale <- scaleallTimeBest

  out <- .GenModusNeu( out$modus, coefs, family$family, X, y, Z, CovRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = 1.05*out$invChol, scale)
  reHess1 <- tcrossprod(out$invChol)

  Us <- .ImportanceSampling(MI = MI2, modus = out$modus, Chol = 1.05*out$invChol, coefs, family$family, X, y, Z, CovRE, w, qvec, ni, scale = scale)

  names(coefs) <- names(beta)

  gr       <- as.numeric( .GradientBetaImp(Us[,-ncol(Us)], coefs, family$family, X, y, Z, w, weightsImp = Us[,ncol(Us)], scale = scale) )
  hess     <- .HessianBetaImp(Us[,-ncol(Us)], coefs, family$family, X, y, Z, w, weightsImp = Us[,ncol(Us)], scale = scale)
  gr2      <- apply(Us[,-ncol(Us)], 1, function(x) tcrossprod(.GradientBetaImp( matrix(x, nrow = 1), coefs, family$family, X, y, Z, w, 1, scale = scale)) )
  gr2      <- colSums( t(gr2) * Us[,ncol(Us)] )
  gr2      <- matrix(gr2, ncol = length(coefs))

  vcovBeta <- hess - gr2 + tcrossprod(gr)
  attr(coefs, "vcov") <- solve( vcovBeta )
  attr(coefs, "vcov.cstr") <- nearPD(solve(vcovBeta))$mat

  # Log Density of the Random Effects Vector
  LLmod <- -.loglikTot(out$modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)
  LLexp <- -.loglikTotImp(Us[,-(n.r+1)], phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, Us[,n.r+1], scale )

  convergence <- max(eps) <= tol1
  if( !convergence ) warning("No convergence!")

  lin.pred <- as.numeric( tcrossprod( X, t(coefs) ) )
  mu <- family$linkinv(lin.pred  + as.numeric(tcrossprod( Z, t(out$modus) ) ))
  names(out$modus) <- colnames(Z)
  REs <- split( out$modus, rep(1:ngrps, ni*qvec) )
  names(REs) <- nam

  resid <- sapply(0:length(REs), function(i) if(i == 0) return(numeric(nrow(X))) else return( as.numeric( crossprod(conLin$reTrms$Ztlist[[i]],  REs[[i]]) ) ) )
  resid <- t( apply(resid, 1, cumsum) )

  REs  <- mapply( function(res, ns){
    m <- matrix(res, nrow = ns, byrow = TRUE)
    rownames(m) <- unique(names(res))
    return(m)}, REs, ni, SIMPLIFY = FALSE )
  nam2 <- names(conLin$reTrms$Ztlist)
  nam2 <- mapply( function(n1, n2) gsub(paste0(" [|] ", n1), "", n2), nam, nam2 )
  nam2 <- gsub("1", "`(Intercept)'", nam2)
  nam2 <- strsplit(nam2, split = " + ", fixed = TRUE)
  REs  <- mapply(function(res, n){
    colnames(res) <- n
    return(res) }, REs, nam2, SIMPLIFY = FALSE)

  attr(REs, "invHess") <- reHess1

  pearson.resid <- (y - mu)/sqrt( family$variance(mu) )
  if(!is.null(weights)) deviance.resid <- family$dev.resids(y, mu, w) else deviance.resid <- family$dev.resids(y, mu, rep(1, length(y)))

  resid <- y - family$linkinv( lin.pred + resid)
  if(!is.null(weights)) sigma <- sqrt( sum( w * resid[,ncol(resid)]^2)/sum(w) ) else sigma <- sqrt( mean( resid[,ncol(resid)]^2) )

  resid <- lapply(1:ncol(resid), function(i) resid[,i])
  attr(resid, "sigma") <- sigma

  V <- CovRE
  names(V) <- nam

  ret <- list(coef = coefs, VarCov = V, scale = scale, RE_mat = Us, RE_mode = REs, residuals = list(pearson = pearson.resid, deviance = deviance.resid, resid = resid), LLmod = LLmod, LLexp = LLexp, niter = i, convergence = convergence, formula = formula, family = family, data = conLin$fr )

  return( ret )

}


