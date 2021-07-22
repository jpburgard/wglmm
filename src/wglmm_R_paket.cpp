// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export(.loglikRE)]]
double loglikRE( List CovRE, arma::vec REs, arma::uvec qvec, arma::uvec ni){

  arma::uword a = qvec.n_elem;
  arma::uword b = 0L;
  arma::uword q;
  arma::uword ns;
  arma::mat V;
  arma::vec us;
  double logDet = 0, sign;
  double logDenRE = 0.0;
  const double logpi = log(M_PI);
  const double log2 = log(2.);

  // Log Density of the Random Effects Vector
  for(uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    V.set_size(q,q);
    V = as<mat>(CovRE[i]);

    log_det(logDet, sign, V);
    logDenRE -= (0.5*logDet*sign*ns);
    logDenRE -= 0.5*q*ns*(log2 + logpi);

    if(det(V) != 0){
      V = inv_sympd(V);
    }else{
      break;
    }

    for(uword j = 0; j < ns; j++) {
      us = REs.subvec(b, b+q-1);
      logDenRE -= 0.5* as_scalar(us.t() * V * us);
      b += q;
    }
  }

  return(logDenRE);
}


// [[Rcpp::export(.loglikCond)]]
double loglikCond(arma::vec REs, arma::vec phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, double scale ){

  arma::uword q;
  arma::uword N = y.n_elem;

  arma::vec eta(N); // Linearer Predictor

  double ret;

  eta = X * phi + Z * REs;

  if(family == "gaussian") {
    double logpi2 = log(2.) + log(M_PI);
    if(weights.n_elem == 1) {
      ret = sum( (eta % y - 0.5*square(eta) - 0.5*square(y)) );
      ret /= scale;
      ret -= y.n_elem * 0.5*( log(scale) + logpi2);
    }else{
      ret = sum( weights % (eta % y - 0.5*square(eta) - 0.5*square(y)) );
      ret /= scale;
      ret -= accu(weights) * 0.5 * ( log(scale) + logpi2);
    }
  }else if(family == "binomial") {
    if(weights.n_elem == 1){
      ret = sum( eta % y - log(1 + exp(eta)) );
    }else{
      ret = sum( weights % (eta % y - log(1 + exp(eta)) ) );
    }
  }else if(family == "poisson") {
    if(weights.n_elem == 1) {
      ret = sum( eta % y - exp(eta) );
    }else{
      ret = sum( weights % (eta % y - exp(eta)) );
    }
  }else if(family == "Gamma" || family == "exponential") {
    if(any(eta <= 0)){
      uvec ind = find(eta <= 0);
      eta( ind ).fill( min(eta(find(eta > 0))) * 0.5 );
    }
    if(weights.n_elem == 1) {
      ret = sum( -scale*eta % y + scale * log(eta) + (scale-1) * log(y) );
      ret += y.n_elem * ( scale * log(scale) - lgamma(scale) );
    }else{
      ret = sum( weights % ( -scale*eta % y + scale * log(eta) + (scale-1) * log(y) ) );
      ret += accu(weights) * ( scale * log(scale) - lgamma(scale) );
    }
  }

  return(ret);
}

// [[Rcpp::export(.loglikTot)]]
double loglikTot( const arma::vec& REs, const arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale){

  double ret;
  ret = loglikCond(REs, phi, family, X, y, Z, weights, scale);
  ret += loglikRE(CovRE, REs, qvec, ni);
  return(-ret);
}

// [[Rcpp::export(.loglikCondImp)]]
double loglikCondImp(const arma::mat& REs, arma::vec phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double scale ){
  arma::uword MI = REs.n_rows;
  double ret = 0;
  for(uword m = 0; m < MI; m++){
    if(weightsImp.at(m) > 0){
      double hv =  loglikCond(REs.row(m).t(), phi, family, X, y, Z, weights, scale);
      ret += weightsImp.at(m) * hv;
    }
  }
  return(ret);
}

// [[Rcpp::export(.loglikREImp)]]
double loglikREImp(List CovRE, const arma::mat& REs, arma::uvec qvec, arma::uvec ni, const arma::vec& weightsImp){
  arma::uword a = qvec.n_elem;
  arma::uword b = 0L;
  arma::uword q;
  arma::uword ns;
  arma::uword MI = REs.n_rows;
  arma::uword m;
  arma::mat V;
  arma::vec us;
  double logDet = 0, sign;
  double logDenRE = 0.0;
  const double logpi = log(M_PI);
  const double log2 = log(2.);

  // Log Density of the Random Effects Vector
  for(uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    V.set_size(q,q);
    V = as<mat>(CovRE[i]);

    log_det(logDet, sign, V);
    logDenRE -= (0.5*logDet*sign*ns);
    logDenRE -= 0.5*q*ns*(log2 + logpi);

    if(det(V) != 0){
      V = inv_sympd(V);
    }else{
      break;
    }

    for(uword j = 0; j < ns; j++) {
      for(m = 0; m < MI; m++){
        us = REs.row(m).subvec(b, b+q-1).t();
        logDenRE -= ( 0.5* as_scalar(us.t() * V * us) * weightsImp.at(m) );
      }
      b += q;
    }
  }
  return(logDenRE);
}

// [[Rcpp::export(.loglikTotImp)]]
double loglikTotImp(const arma::mat& REs, const arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, const arma::vec& weightsImp, double scale){

  double ret;

  ret = loglikCondImp(REs, phi, family, X, y, Z, weights, weightsImp, scale);
  ret += loglikREImp(CovRE, REs, qvec, ni, weightsImp);

  return(-ret);
}


arma::vec Gradient(int b, int q, int k, arma::vec &preCalcValues, const arma::vec& REs, const List& invmats){
  arma::vec Grad(q, fill::zeros);

  for(uword i = 0; i < q; i++) {
    Grad.at(i) = preCalcValues.at(b+i);
  }

  arma::mat M = as<arma::mat>(invmats[k]);

  Grad = Grad - M * REs.subvec(b,b+q-1);

  return(Grad);
}

// [[Rcpp::export(.GradientAll)]]
arma::vec GradientAll(const arma::vec& REs, const arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale){
  arma::uword nr = REs.n_elem;
  arma::uword a = qvec.n_elem;
  arma::uword b = 0;
  arma::uword q;
  arma::uword ns;
  arma::vec grad(nr);

  arma::vec eta = X * phi + Z * REs;

  List invmats(a);
  arma::mat M;
  for(q = 0; q < a; q++){
    M = as<mat>(CovRE[q]);
    M = inv_sympd(M);
    invmats(q) = M;
  }

  // ------------------------------------- Baue Gradienten ---------------------------------------------------------
  if(family == "binomial") {
    eta = exp(eta);
    eta = eta/(1+eta);
  }else if(family == "poisson") {
    eta = exp(eta);
  }else if(family == "Gamma" || family == "exponential") {
    if(any(eta <= 0)){
      uvec ind = find(eta <= 0);
      eta( ind ).fill( min(eta(find(eta > 0))) * 0.5 );
    }
    eta = 1/eta;
  }

  //Vorberechnung für die Werte für Gradient. Muss nur einmal gemacht werden.
  arma::vec yMinuseta = y - eta;
  if(family == "gaussian") yMinuseta /= scale;
  if(family == "Gamma") yMinuseta *= -scale;

  if(weights.n_elem > 1){
    yMinuseta = weights % yMinuseta;
  }

  arma::vec preCalcValues(Z.n_cols, fill::zeros);

  for(uword i = 0; i < Z.n_cols; i++) {
    preCalcValues.at(i) = accu(yMinuseta % Z.col(i));
  }

  for( uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    for(uword j = 0; j < ns; j++) {
      grad.subvec(b, b+q-1) = Gradient(b, q, i, preCalcValues, REs, invmats);
      b +=q;
    }
  }

  return(-grad);
}

// [[Rcpp::export(.GradientBetaImp)]]
arma::vec GradientBetaImp(const arma::mat& REs, const arma::vec& phi, const String& family, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double scale){

  arma::uword N = y.n_elem;
  arma::uword ngrad = phi.n_elem;
  arma::uword MI = weightsImp.n_elem;

  arma::vec eta(N); // Linear Predictor
  arma::vec Grad(ngrad, fill::zeros);

  int XRows = X.n_rows;
  arma::vec Xphi = X * phi;
  bool weightsIn = weights.n_elem > 1;

  for(uword m = 0; m < MI; m++) {

    double weightsM = weightsImp.at(m);

    eta = Xphi + Z * REs.row(m).t();

    //  --------------------------------------- Baue Gradienten ------------

    if(family == "binomial") {
      eta = exp(eta);
      eta = eta/(1+eta);
    }else if(family == "poisson") {
      eta = exp(eta);
    }else if(family == "Gamma" || family == "exponential") {
      if(any(eta <= 0)){
        uvec ind = find(eta <= 0);
        eta( ind ).fill( min(eta(find(eta > 0))) * 0.5 );
      }
      eta = 1/eta;
    }
    arma::vec yMinuseta = y - eta;
    if(family == "gaussian") yMinuseta /= scale;
    if(family == "Gamma" || family == "exponential") {
      yMinuseta *= -scale;
    }
    if(weightsIn){
      yMinuseta = weights % yMinuseta;
    }

    for(uword i = 0; i < ngrad; i++) {
      double* mem = X.colptr(i);
      arma::vec Xcol(mem, XRows, false, false);
      Grad.at(i) += weightsM * accu( yMinuseta % Xcol );
    }

  }

  return(-Grad); // negative ll is evaluated

}

// [[Rcpp::export(.HessianBetaImp)]]
arma::mat HessianBetaImp(const arma::mat& REs, const arma::vec& phi, const String family, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double scale){

  arma::uword nr = Z.n_rows;
  arma::uword N = y.n_elem;
  arma::uword nhess = phi.n_elem;
  arma::uword MI = REs.n_rows;
  arma::vec hv(N);
  arma::mat Hess(nhess, nhess, fill::zeros);

  bool weightsIn = weights.n_elem > 1;

  if(family != "gaussian"){
    arma::vec Xphi = X * phi;
    arma::vec eta(N);

    for(uword m = 0; m < MI; m++) {

      // -------------------------------------- Linear predictors with REs ------------------
      eta = Xphi + Z * REs.row(m).t();
      //  --------------------------------------- Hessian --------------------------------------------------
      if(family == "binomial") {
        eta = exp(eta);
        eta = eta/square(1+eta);
      }else if(family == "poisson") {
        eta = -exp(eta);
      }else if(family == "Gamma" || family == "exponential"){
        if(any(eta <= 0)){
          uvec ind = find(eta <= 0);
          eta( ind ).fill( min(eta(find(eta > 0))) * 0.5 );
        }
        eta = scale/square(eta);
      }

      if(weightsIn){
        eta = eta % weights;
      }

      double* mem;
      double* mem2;
      double weightsM = weightsImp.at(m);
      int rows = X.n_rows;

      for(uword i = 0; i < nhess; i++) {
        mem = X.colptr(i);
        arma::vec Xcoli(mem, rows, false, false);
        hv = eta % Xcoli;
        for(uword j = 0; j <= i; j++) {
          mem2 = X.colptr(j);
          arma::vec Xcolj(mem2, rows, false, false);
          Hess.at(i,j) -= weightsM*accu( hv % Xcolj );
        }
      }
    }
  }else{
    double* mem;
    double* mem2;
    int rows = X.n_rows;

    for(uword i = 0; i < nhess; i++) {
      mem = X.colptr(i);
      arma::vec Xcoli(mem, rows, false, false);
      if(weightsIn){
        hv = Xcoli % weights;
      }else{
        hv = Xcoli;
      }
      for(uword j = 0; j <= i; j++) {
        mem2 = X.colptr(j);
        arma::vec Xcolj(mem2, rows, false, false);
        Hess.at(i,j) -= accu( hv % Xcolj )/scale;
      }
    }

  }


  Hess = arma::symmatl(Hess);

  return(-Hess);

}



// // #####################################################################################################################
// // ################# Newton Raphson to Find New Coefficients ###########################################################
// // #####################################################################################################################

// [[Rcpp::export(.NewCoefImp)]]
List NewCoefImp(const arma::mat& REs, arma::vec& phi, const String& family, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double tol, int MaxIt, double scale ){

  arma::uword ncoef = phi.n_elem;
  arma::uword ind;
  arma::uvec ind2;
  arma::vec phiNew(ncoef);
  arma::vec phiOld = phi;
  arma::mat Hess(ncoef, ncoef);
  arma::vec GradOld(ncoef);
  arma::vec dir(ncoef);

  double eps = 1.1 * tol;

  bool klappt;
  List ret;

  for(int a = 0; a < MaxIt; a++) {

    Hess = HessianBetaImp(REs, phiOld, family, X, y, Z, weights, weightsImp, scale);
    GradOld = GradientBetaImp(REs, phiOld, family, X, y, Z, weights, weightsImp, scale);
    klappt = solve(dir, Hess, GradOld); //, solve_opts::no_approx
    if(klappt) {
      phiNew = phiOld - dir;
    }else{
      break;
    }

    eps = max( abs( phiNew - phiOld)/( abs(phiOld) + 1) );
    phiOld = phiNew;
    if(eps < tol) {
      break;
    }
  }

  if(eps > tol) {
    Rcpp::Rcout << "Warning: No convergence of NR!" << endl;
    klappt = false;
  }

  ret("phi") = phiNew;
  ret("convergence") = klappt;

  return(ret);
}


// [[Rcpp::export(.NewCovREImp)]]
List NewCovREImp( const arma::mat& REs, const arma::uvec& qvec, const arma::uvec& ni, const arma::vec& weightsImp){
  arma::uword MI = weightsImp.n_rows;
  arma::uword ns;
  arma::uword q;
  arma::uword b = 0;
  int a = qvec.n_elem;
  arma::mat N;
  arma::mat us;
  bool tries;

  List newCovRE(a);

  for(uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    N.set_size(q,q);
    N.zeros();

    us.set_size(q, MI);

    double* mem;
    int usRows = q;

    for(uword j = 0; j < ns; j++) {
      us = REs.cols(b, b+q-1).t();

      for(uword m = 0; m < MI; m++) {
        mem = us.colptr(m);
        vec usCol(mem, usRows, false, true);
        N += (usCol * usCol.t() * as_scalar( weightsImp.at(m)) );
      }
      b += q;
    }

    N /= ns;
    newCovRE[i] = N;
  }

  return(newCovRE);
}


// //#######################################################################################################################
// //#######################################################################################################################
// //######################### Samplers ####################################################################################
double evalhv( const arma::vec& REs, const List& invmat, const arma::uvec& qvec, const arma::uvec& ni){
  arma::uword a = qvec.n_elem;
  arma::uword b = 0;
  arma::uword q, ns;
  arma::uword i, j;
  double ret = 0.;

  for(i = 0; i < a; i++){
    q = qvec.at(i);
    ns = ni.at(i);
    arma::mat M = as<arma::mat>(invmat[i]);
    vec us(q);

    for(j = 0; j < ns; j++){
      us = REs.subvec(b, b+q-1);
      ret -= 0.5*as_scalar( us.t()*M*us );
      b += q;
    }
  }
  return(ret);
}

// [[Rcpp::export(.ImportanceSampling)]]
arma::mat ImportanceSampling( int MI, const arma::vec& modus, const arma::mat& Chol, arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale){

  arma::uword nr = Z.n_cols;
  arma::mat rands(nr, MI);
  List control;
  arma::vec re(nr);
  arma::vec us(nr);
  arma::vec w(MI, fill::zeros);
  double sup;
  double sign;
  double ldet;
  double addLdet;

  sup = loglikCond( modus, phi, family, X, y, Z, weights, scale);

  arma::uword q;
  arma::uword a = qvec.n_elem;
  double hv;
  List invmat(a);

  for(uword i = 0; i < a; i++){
    q = qvec.at(i);
    arma::mat M = as<arma::mat>(CovRE[i]);
    if( q > 1){
      invmat(i) = arma::inv_sympd(M);
    }else{
      invmat(i) = 1/M;
    }
  }

  for(uword k = 0; k < MI; k = k+2) {
    us.randn();
    addLdet = 0.5*as_scalar( us.t() * us ); // If V inverse of covariance matrix and L*L.t() = V^{-1}, then re.t() * L.t() * V * L * re = re.t() * re
    re = Chol * us;
    arma::vec randsCol = modus + re;
    rands.col(k) = randsCol;
    hv = evalhv( randsCol, invmat, qvec, ni);
    w.at(k) += addLdet;
    w.at(k) += loglikCond(randsCol, phi, family, X, y, Z, weights, scale);
    w.at(k) += hv;
    randsCol = modus - re;
    rands.col(k+1) = randsCol;
    hv = evalhv( randsCol, invmat, qvec, ni);
    w.at(k+1) += addLdet;
    w.at(k+1) += loglikCond(randsCol, phi, family, X, y, Z, weights, scale);
    w.at(k+1) += hv;
  }

  rands = rands.t();
  uvec ind = find( w == w );
  w = w(ind);
  rands = rands.rows(ind);
  w -= sup;
  w -= mean(w);
  int count = 0;

  for( uword i = 0; i < w.n_elem; i++) {
    double value = w.at(i);
    if (value > 35) {
      // cout << "A" << endl;
      w.at(i) = 35;
      count++;
    } else if (value < -25) {
      // cout << "B" << endl;
      w.at(i) = -25;
      count++;
    }
  }

  w = arma::exp(w);
  w /= arma::sum(w);

  rands = arma::join_horiz(rands, w);

  return(rands);

}


// BFGS for modes of Random Effects
// [[Rcpp::export(.GenModusNeu)]]
List GenModusNeu(const arma::vec& REs, const arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double tol, int MaxIt, arma::mat Chol, double scale){

  arma::uword nr = REs.n_elem;
  arma::uword ind;
  arma::uvec ind2;
  arma::vec modusNeu(nr);
  arma::vec modusAlt(nr);
  arma::vec gradAlt(nr);
  arma::vec gradNeu(nr);
  arma::vec dir(nr);
  arma::vec res(nr);
  arma::vec alpha = linspace<vec>(0.001,1,45);
  arma::vec ll(45);
  arma::vec h(nr);
  arma::vec s(nr);

  double eps;
  double factor;
  double phifac;
  double sigmaNeu;
  double sigmaAlt = 1.;
  bool rescale;
  bool conv = false;
  double LLMin;
  int k;
  List ret;

  modusAlt = REs;
  gradAlt = GradientAll(modusAlt, phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);
  dir = -Chol * Chol.t() * gradAlt;
  LLMin = loglikTot(modusAlt, phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);

  double* mem;
  int CholRows = Chol.n_rows;

  for(k = 0; k < MaxIt; k++) {
    for(uword i = 0; i < 45; i++){
      double hv = loglikTot(modusAlt + dir*alpha.at(i), phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);
      if(hv == hv){
        ll.at(i) = hv;
      }else{
        ll.at(i) == LLMin;
      }
    }
    if(LLMin <= min(ll)){
      alpha *= 0.1;
      for(uword i = 0; i < 45; i++){
        double hv = loglikTot(modusAlt + dir*alpha.at(i), phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);
        if(hv == hv){
          ll.at(i) = hv;
        }else{
          ll.at(i) == LLMin;
        }
      }
    }

    if(LLMin <= min(ll)){
      alpha *= 0.1;
      for(uword i = 0; i < 45; i++){
        double hv = loglikTot(modusAlt + dir*alpha.at(i), phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);
        if(hv == hv){
          ll.at(i) = hv;
        }else{
          ll.at(i) == LLMin;
        }
      }
    }

    if(LLMin <= min(ll)) break;
    LLMin = arma::min(ll);
    ind2 = arma::sort_index(ll);
    ind = ind2(0);
    dir = dir* alpha.at(ind);
    modusNeu = modusAlt + dir;
    gradNeu = GradientAll(modusNeu, phi, family, X, y, Z, CovRE, weights, qvec, ni, scale);
    res = gradNeu - gradAlt;

    if( norm(modusNeu-modusAlt) < tol ) {
      conv = true;
      break;
    }
    factor = 1/accu( res % dir);

    s = Chol.t() * res;
    ind2 = find(s, 1, "last");
    ind = ind2(0L);

    modusAlt = modusNeu;

    h = s.at(ind) * Chol.col(ind);
    phifac = s.at(ind)*s.at(ind);

    for(uword i = ind; i > 0; i--) {
      mem = Chol.colptr(i-1);
      arma::vec CholCol(mem, CholRows, false, false);
      vec hvvec = sqrt(phifac/(s.at(i-1)*s.at(i-1) + phifac)) * ( s.at(i-1)/phifac * h - CholCol);
      if(hvvec.has_nan()) break;
      Chol.col(i) = hvvec;

      phifac += s(i-1)*s.at(i-1);
      h += s.at(i-1)*CholCol;
    }

    Chol.col(0L) = dir * sqrt(factor);

    sigmaNeu = min( sigmaAlt, norm(Chol.col(0L)) );
    rescale = sigmaNeu == sigmaAlt;

    Chol.cols(1L, nr-1) = Chol.cols(1, nr-1) - dir * (res.t() * Chol.cols(1,nr-1) * factor);
    if(rescale) {
      Chol.cols(1, nr-1) = normalise(Chol.cols(1L, nr-1L)) * sigmaNeu;
    }

    sigmaAlt = sigmaNeu;
    // approx. Cholesky of INVERSE Hessian
    modusAlt = modusNeu;
    gradAlt = gradNeu;
    dir = -Chol*Chol.t()*gradAlt;
  }

  if(!conv) {
    Rcpp::Rcout << "No convergence!" << endl;
  }
  ret("modus") = modusNeu;
  ret("invChol") = Chol;
  ret("iter") = k;

  return(ret);
}

double likImp(int MI, const arma::vec& phi, const String& family, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale){
  arma::uword a = qvec.n_elem;
  arma::uword nr = Z.n_cols;
  arma::uword b = 0;
  arma::uword i, j, m, q, ns;
  arma::vec us;
  arma::mat REs(nr,MI);
  arma::mat M;
  double likeli = 0.;

  double scaling = loglikCond(zeros(nr), phi, family, X, y, Z, weights, scale);

  for( i = 0; i < a; i++){
    q = qvec.at(i);
    ns = ni.at(i);
    us.set_size(q);
    M = as<arma::mat>(CovRE[i]);
    M = chol(M);

    for(j = 0; j < ns; j++){
      for(m = 0; m < MI; m=m+2){
        us.randn();
        us = M * us;
        REs.col(m).subvec(b,b+q-1) = us;
        REs.col(m+1).subvec(b,b+q-1) = -us;
      }
      b += q;
    }
  }

  us.set_size(nr);

  for( m = 0; m < MI; m++){
    us = REs.col(m);
    double ll = loglikCond(us, phi, family, X, y, Z, weights, scale);
    likeli += exp(ll - scaling);
  }
  likeli /= MI;
  likeli = log(likeli) + scaling;

  return(likeli);
}

