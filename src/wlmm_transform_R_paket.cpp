// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

arma::vec trans(const arma::vec& y, const String& transform, double lambda){
  arma::vec ytilde = y;
  if( (transform != "none") & (lambda != 0) ){
    ytilde = arma::pow(y, lambda );
    if(transform == "dual"){
      ytilde = 0.5*(ytilde - 1/ytilde)/lambda;
    }else if(transform == "box-cox"){
      ytilde--;
      ytilde /= lambda;
    }
  }else if((transform != "none") & (lambda == 0)){
    ytilde = log(ytilde);
  }
  return(ytilde);
}


arma::vec jacobian_transform(const arma::vec& y, const String& transform, double lambda){
  arma::vec ytilde = y;
  if(transform != "none"){
    ytilde = arma::pow(y, lambda );
    if(transform == "dual"){
      ytilde = 0.5*(ytilde + 1/ytilde)/y;
    }else if(transform == "box-cox"){
      ytilde /= y;
    }
  }
  return(ytilde);
}

double loglikRE_trans( List CovRE, arma::vec REs, arma::uvec qvec, arma::uvec ni){

  arma::uword a = qvec.n_elem;
  arma::uword b = 0L;
  arma::uword q;
  arma::uword ns;
  arma::mat V;
  arma::vec us;
  double logDet = 0, sign;
  double logDenRE = 0.0;

  // Log Density of the Random Effects Vector
  for(arma::uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    V.set_size(q,q);
    V = as<arma::mat>(CovRE[i]);

    arma::log_det(logDet, sign, V);
    logDenRE -= (0.5*logDet*sign * ns);

    if(q == 1){
      V = 1/V;
    }else if(!(q == 1 && arma::all(V.diag() == 0))) {
      V = arma::inv_sympd(V);
    }else{
      V = arma::pinv(V);
    }

    for(arma::uword j = 0; j < ns; j++) {
      us = REs.subvec(b, b+q-1);
      logDenRE -= 0.5* arma::as_scalar(us.t() * V * us);
      b += q;
    }
  }

  return(logDenRE);
}


double loglikCond_trans( const arma::vec& REs, arma::vec phi, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, double scale ){ //, const uvec& qvec, const uvec& ni) {

  arma::uword q;
  arma::uword N = y.n_elem;

  bool weightsIn = (weights.n_elem > 1);

  arma::vec eta(N); // Linearer Praediktor

  double ret;

  eta = X * phi + Z * REs;

  if(weightsIn) {
    ret = sum( weights % (eta % y - 0.5*square(eta) - 0.5*square(y)) );
    ret /= scale;
    ret -= accu(weights) * 0.5*log(scale);
  }else{
    ret = sum( (eta % y - 0.5*square(eta) - 0.5*square(y)) );
    ret /= scale;
    ret -= N * 0.5*log(scale);
  }

  return(ret);
}

double loglik_lambda( const arma::vec& linpredz, const String& transform, const arma::mat& X, const arma::vec& y, const arma::vec& weights, const double& lambda, const arma::mat& xwxxw){

  double sigmahat;
  bool weightsIn = weights.n_elem > 1;
  arma::uword N = y.n_elem;
  arma::vec ytrans = trans(y, transform, lambda);
  double ll = 0.;

  arma::vec betahat = xwxxw * (ytrans - linpredz); // solve_opts::likely_sympd
  arma::vec linpred = X * betahat + linpredz;
  arma::vec res = ytrans - linpred;
  if(weightsIn){
    double sumw = accu(weights);
    sigmahat = accu(weights % square(res))/sumw;
  }else{
    sigmahat = accu( square(res))/N;
  }

  ll -= 0.5 * N * log(sigmahat);
  if(weightsIn){
    ll += arma::accu(weights % arma::log(jacobian_transform(y, transform, lambda)));
  }else{
    ll += arma::accu( arma::log(jacobian_transform(y, transform, lambda)));
  }

  return(-ll);
}

double logliktTot_trans( const arma::vec& REs, const arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale, double lambda, const mat& xwxxw){

  double ret;
  arma::vec ytrans = trans(y, transform, lambda);
  // ret = -loglik_lambda( Z * REs, transform, X, y, weights, lambda, xwxxw);
  ret = loglikCond_trans( REs, phi, X, ytrans, Z, weights, scale);
  ret += loglikRE_trans(CovRE, REs, qvec, ni);
  if(weights.n_elem == 1){
    ret += accu(log(jacobian_transform(y, transform, lambda)));
  }else{
    ret += accu(weights % log(jacobian_transform(y, transform, lambda)));
  }
  return(-ret);
}

double loglikCondImp_trans( const arma::mat& REs, arma::vec phi, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double scale ){
  arma::uword MI = REs.n_rows;
  double ret = 0;
  for(arma::uword m = 0; m < MI; m++){
    arma::vec hv = REs.row(m).t();
    ret += weightsImp.at(m) * loglikCond_trans( hv, phi, X, y, Z, weights, scale);
  }
  return(ret);
}

double loglikREImp_trans(List CovRE, const arma::mat& REs, arma::uvec qvec, arma::uvec ni, const arma::vec& weightsImp){
  arma::uword MI = REs.n_cols;
  arma::uword a = qvec.n_elem;
  arma::uword b = 0L;
  arma::uword q;
  arma::uword ns;
  arma::mat V;
  arma::vec us;
  double logDet = 0, sign;
  double logDenRE = 0.0;

  for(arma::uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    V.set_size(q,q);
    V = as<mat>(CovRE[i]);

    log_det(logDet, sign, V);
    logDenRE -= (0.5*logDet*sign *ns);

    if(!(q == 1 && all(V.diag() == 0))) {
      V = arma::inv_sympd(V);
    }else{
      V = pinv(V);
    }
  }
  for(arma::uword j = 0; j < ns; j++) {
    for(arma::uword m = 0; m < MI; m++){
      us = REs.row(m).subvec(b, b+q-1).t();
      logDenRE -= 0.5* as_scalar(us.t() * V * us) * weightsImp(m);
    }
    b += q;
  }

  return(logDenRE);
}

double loglik_lambdaImp( const arma::mat& REs, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat Z, const arma::vec& weights, const arma::vec& weightsImp, const double& lambda, const arma::mat& xwxxw){

  double ll = 0.;

  arma::uword N = y.n_elem;
  arma::uword m, MI = REs.n_rows;
  double sigmahat = 0;
  arma::vec betahat(X.n_cols);
  arma::vec ytrans = trans(y, transform, lambda);
  arma::vec linpredz(N);

  bool weightsIn = (weights.n_elem > 1);

  linpredz = Z * sum( REs.each_col() % weightsImp, 0).t();
  betahat = xwxxw * (ytrans - linpredz);
  arma::vec linpred = X * betahat;
  arma::vec res = ytrans - linpred;

  if(weightsIn){
    for(m = 0; m < MI; m++){
      arma::vec us = REs.row(m).t();
      sigmahat += weightsImp.at(m) * accu( weights % square(res - Z * us));
    }
    // sigmahat = accu(weights % square(res - linpredz));
  }else{
    for(m = 0; m < MI; m++){
      arma::vec us = REs.row(m).t();
      sigmahat += weightsImp.at(m) * accu( square(res - Z * us));
    }
    // sigmahat = accu( square(res - linpredz));
  }
  sigmahat /= N;

  ll -= 0.5 * N * log(sigmahat);
  if(weightsIn){
    ll += arma::accu(weights % arma::log(jacobian_transform(y, transform, lambda)));
  }else{
    ll += arma::accu( arma::log(jacobian_transform(y, transform, lambda)));
  }

  return(-ll);
}


double logliktTotImp_trans(const arma::mat& REs, const arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, const arma::vec& weightsImp, double scale, double lambda, const arma::mat& xwxxw){
  uword MI = REs.n_rows;
  uword ns;
  uword q;
  uword a = qvec.n_elem;
  double ret = 0;

  arma::vec ytrans = trans(y, transform, lambda);
  arma::vec linpredz(y.n_elem, fill::zeros);
  // ret = -loglik_lambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambda, xwxxw);
  ret = loglikCondImp_trans( REs, phi, X, ytrans, Z,  weights, weightsImp, scale);
  ret += loglikREImp_trans(CovRE, REs, qvec, ni, weightsImp);

  if(weights.n_elem == 1){
    ret += accu(log(jacobian_transform(y, transform, lambda)));
  }else{
    ret += accu(weights % log(jacobian_transform(y, transform, lambda)));
  }
  return(-ret);
}

vec Gradient_trans(int b, int q, int k, arma::vec &preCalcValues, const arma::vec& REs, const List& invmats){
  arma::vec Grad(q, fill::zeros);

  for(arma::uword i = 0; i < q; i++) {
    Grad.at(i) = preCalcValues.at(b+i);
  }

  arma::mat M = as<arma::mat>(invmats[k]);

  Grad = Grad - M * REs.subvec(b,b+q-1);

  return(Grad);
}

vec GradientAll_trans(const arma::vec& REs, const arma::vec& phi, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale, const arma::mat& xwxxw){
  arma::uword nr = REs.n_elem;
  arma::uword a = qvec.n_elem;
  arma::uword b = 0;
  arma::uword q;
  arma::uword ns;

  List invmats(a);
  arma::mat M;
  for(q = 0; q < a; q++){
    M = as<arma::mat>(CovRE[q]);
    if(M.n_rows > 1) {
      M= arma::inv_sympd(M);
    }else{
      if((as_scalar(M) != 0)) M = 1/M;
    }
    invmats(q) = M;
  }

  bool weightsIn = (weights.n_elem > 1);
  arma::vec beta = phi;
  double newscale = scale;

  // // special
  // beta = xwxxw * (y - Z * REs);
  // if(weightsIn){
  //   newscale = accu(weights % square(y - X * beta - Z * REs))/y.n_elem;
  // }else{
  //   newscale = accu(square(y - X * beta - Z * REs))/y.n_elem;
  // }

  arma::vec grad(nr);

  arma::vec eta = X * beta + Z * REs;

  // ------------------------------------- Baue Gradienten --------------------------------------------------------

  //Vorberechnung für die Werte für Gradient. Muss nur einmal gemacht werden.
  arma::vec yMinuseta = y - eta;
  yMinuseta /= newscale;

  if(weightsIn){
    yMinuseta = weights % yMinuseta;
  }

  arma::vec preCalcValues(Z.n_cols, fill::zeros);

  for(arma::uword i = 0; i < Z.n_cols; i++) {
    // preCalcValues.at(i) = accu(yMinuseta % (Z.col(i) - X * xwxxw * Z.col(i)) );
    preCalcValues.at(i) = accu(yMinuseta % Z.col(i));
  }

  for( arma::uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    for(arma::uword j = 0; j < ns; j++) {
      grad.subvec(b, b+q-1) = Gradient_trans(b, q, i, preCalcValues, REs, invmats);
      b +=q;
    }
  }

  return(-grad);
}

vec GradientBetaImp_trans( const arma::mat& REs, const arma::vec& phi, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double scale){

  arma::uword N = y.n_elem;
  arma::uword ngrad = phi.n_elem;
  arma::uword MI = weightsImp.n_elem;

  arma::vec eta(N); // Linearer Praediktor
  arma::vec Grad(ngrad, fill::zeros);

  int XRows = X.n_rows;
  arma::vec Xphi = X * phi;

  for(arma::uword m = 0; m < MI; m++) {

    double weightsM = weightsImp.at(m);

    eta = Xphi + Z * REs.row(m).t();

    //  --------------------------------------- Baue Gradienten --------------------------------------------------
    arma::vec yMinuseta = y - eta;

    if(weights.n_elem > 1){
      yMinuseta = weights % yMinuseta;
    }

    for(arma::uword i = 0; i < ngrad; i++) {
      double* mem = X.colptr(i);
      arma::vec Xcol(mem, XRows, false, false);
      Grad.at(i) += weightsM * arma::sum( yMinuseta % Xcol )/scale;
    }

  }

  return(-Grad); // Es soll ja die negative LL evaluiert werden

}

arma::mat HessianBetaImp_trans( arma::mat& X, const arma::vec& weights, double scale){

  arma::uword N = X.n_rows;
  arma::uword nhess = X.n_cols;
  bool weightsIn = weights.n_elem > 1;
  arma::vec eta(N);
  arma::mat Hess(nhess, nhess, fill::zeros);

  int rows = X.n_rows;


  for(arma::uword i = 0; i < nhess; i++) {
    double* mem;
    mem = X.colptr(i);
    arma::vec Xcoli(mem, rows, false, false);
    for(arma::uword j = 0; j <= i; j++) {
      double* mem2;
      mem2 = X.colptr(j);
      arma::vec Xcolj(mem2, rows, false, false);
      if(weightsIn){
        Hess.at(i,j) -= sum( Xcoli % Xcolj % weights ) / scale;
      }else{
        Hess.at(i,j) -= sum( Xcoli % Xcolj ) / scale;
      }

    }
  }

  Hess = symmatl(Hess);

  return(Hess);

}

double GradientLambdaImp( const arma::mat& REs, const String& transform, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, const double lambda, const arma::mat& xwxxw){

  double Grad = 0.;
  double Gradhv = 0.;

  arma::uword N = y.n_elem;
  arma::vec ytrans = trans(y, transform, lambda);
  arma::vec trans_lambda1(N);
  arma::vec logy = log(y);
  bool weightsIn = weights.n_elem > 1;

  arma::vec linpredz = Z * sum( REs.each_col() % weightsImp, 0).t();
  arma::vec phi = xwxxw * (ytrans - linpredz);
  arma::vec linpred = X * phi;
  arma::vec res = ytrans - linpred;

  double scale = 0;
  if(weightsIn){
    double sumw = accu(weights);
    for(arma::uword m = 0; m < REs.n_rows; m++){
      arma::vec us = REs.row(m).t();
      scale += weightsImp.at(m) * accu( weights % square(res - Z * us) )/sumw;
    }
  }else{
    for(arma::uword m = 0; m < REs.n_rows; m++){
      arma::vec us = REs.row(m).t();
      scale += weightsImp.at(m) * accu( square(res - Z * us) )/N;
    }
  }

  // arma::vec us = sum(REs.each_col() % weightsImp, 0).t();
  // if(weightsIn){
  //     scale += accu( weights % square(res - Z * us) )/N;
  // }else{
  //   scale += accu( square(res - Z * us) )/N;
  // }

  if( transform == "box-cox"){
    if(lambda == 0){
      trans_lambda1 = square(logy)*0.5;
    }else{
      trans_lambda1 = (pow(y, lambda) % logy - ytrans)/lambda;
    }
    if(weightsIn){
      Grad = accu(weights % logy);
    }else{
      Grad = accu(logy);
    }
  }else if(transform == "dual"){
    arma::vec ytrans_jacob = jacobian_transform(y, transform, lambda);
    if(lambda == 0){
      trans_lambda1.zeros();
    }else{
      trans_lambda1 = (logy % ytrans_jacob % y - ytrans)/lambda;
    }
    if(weightsIn){
      Grad = accu( weights % logy % ytrans / (y % ytrans_jacob) ) * lambda;
    }else{
      Grad = accu( logy %  ytrans / ( y % ytrans_jacob ) ) * lambda;
    }
  }

  res -= linpredz;
  // if(weightsIn){
  //   Grad -= accu( weights % trans_lambda1 % res) /scale;
  // }else{
  //   Grad -= accu( trans_lambda1 % res) /scale;
  // }
  arma::vec resderiv = trans_lambda1 - X * xwxxw * trans_lambda1;
  if(weightsIn){
    Grad -= arma::accu( weights % resderiv % res) /scale;
  }else{
    Grad -= arma::accu( resderiv % res) /scale;
  }

  return(-Grad); // Negative LL
}

double HessianLambdaImp( const arma::mat& REs, const String& transform, const arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double lambda, const arma::mat& xwxxw){
  arma::uword N = y.n_elem;

  arma::vec ytrans = trans(y, transform, lambda);
  double Hess = 0.;
  arma::vec trans_lambda2(N);
  arma::vec trans_lambda1(N);
  arma::vec logy = log(y);
  bool weightsIn = weights.n_elem > 1;

  arma::vec linpredz = Z * sum( REs.each_col() % weightsImp, 0).t();
  arma::vec phi = xwxxw * (ytrans - linpredz);
  arma::vec linpred = X * phi;
  arma::vec res = ytrans - linpred;
  double scale = 0;
  if(weightsIn){
    for(arma::uword m = 0; m < REs.n_rows; m++){
      arma::vec us = REs.row(m).t();
      scale += weightsImp.at(m) * accu( weights % square(res - Z * us) );
    }
  }else{
    for(arma::uword m = 0; m < REs.n_rows; m++){
      arma::vec us = REs.row(m).t();
      scale += weightsImp.at(m) * accu( square(res - Z * us) );
    }
  }
  scale /= N;

  // arma::vec us = sum(REs.each_col() % weightsImp, 0).t();
  // if(weightsIn){
  //   scale += accu( weights % square(res - Z * us) )/N;
  // }else{
  //   scale += accu( square(res - Z * us) )/N;
  // }

  res = res - linpredz;

  if( transform == "box-cox"){
    if(lambda != 0){
      arma::vec ylam = pow(y, lambda);
      trans_lambda1 = (logy % ylam - ytrans)/lambda;
      trans_lambda2 = square(logy) % ylam / lambda - 2 * logy % ylam / (lambda * lambda) + 2 * ytrans/(lambda * lambda);
    }else if(lambda == 0){
      trans_lambda1 = 0.5 * square(logy);
      trans_lambda2 = (2. / 3) * trans_lambda1 % logy;
    }

  }else if(transform == "dual"){
    arma::vec ytrans_jacob = jacobian_transform(y, transform, lambda);
    if(lambda != 0 ){
      trans_lambda1 = (logy % ytrans_jacob % y - ytrans) / lambda;
      trans_lambda2 = square(logy) % ytrans - 2 * logy % ytrans_jacob % y/(lambda * lambda) + 2*ytrans/(lambda*lambda);
    }else if(lambda == 0){
      trans_lambda1.zeros();
      trans_lambda2 = square( logy) % logy / 3;
    }

    if(weightsIn){
      Hess += arma::accu(weights %  arma::square(logy) % (1 - lambda * lambda * arma::square( ytrans/( y % ytrans_jacob) ) ) );
    }else{
      Hess += arma::accu( arma::square(logy) % ( 1 - lambda * lambda * arma::square( ytrans/(y % ytrans_jacob) ) ) );
    }
  }

  // arma::vec prefac1 = (res % trans_lambda2) + square(trans_lambda1);
  //
  // if(weightsIn){
  //   double prefac2 = accu(weights % res % trans_lambda1);
  //   prefac2 *= prefac2 * 2.;
  //   prefac2 /= (N * scale * scale);
  //   Hess -= accu(weights % prefac1)/scale;
  //   Hess += prefac2;
  // }else{
  //   double prefac2 = accu(res % trans_lambda1);
  //   prefac2 *= prefac2 * 2.;
  //   prefac2 /= (N * scale * scale);
  //   Hess -= accu(prefac1)/scale;
  //   Hess += prefac2;
  // }

  arma::vec resderiv = trans_lambda1 - X * (xwxxw * trans_lambda1);
  arma::vec prefac1 = (res % trans_lambda2) + (resderiv % trans_lambda1);

  if(weightsIn){
    double prefac2 = accu(weights % res % resderiv);
    prefac2 *= accu(weights % res % trans_lambda1) * 2.;
    prefac2 /= (N * scale * scale);
    Hess -= accu(weights % prefac1)/scale;
    Hess += prefac2;
  }else{
    double prefac2 = accu(res % resderiv);
    prefac2 *= arma::accu(res % trans_lambda1) * 2.;
    prefac2 /= (N * scale * scale);
    Hess -= accu(prefac1)/scale;
    Hess += prefac2;
  }

  // arma::vec resderiv = trans_lambda1 - X * (xwxxw * trans_lambda1);
  // arma::vec resderiv2 = trans_lambda2 - X * xwxxw * trans_lambda2;
  // if(weightsIn){
  //   double hv = accu( weights % res % resderiv)/scale;
  //   Hess += hv * hv;
  //   Hess -= accu( weights % square(resderiv))/scale;
  //   Hess -= accu(weights % res % resderiv2)/scale;
  // }else{
  //   double hv = accu( res % resderiv)/scale;
  //   Hess += hv * hv;
  //   Hess -= accu( square(resderiv))/scale;
  //   Hess -= accu(res % resderiv2)/scale;
  // }

  return(-Hess);
}

// // #####################################################################################################################
// // ################# Newton Raphson to Find New Coefficients ###########################################################
// // #####################################################################################################################
arma::vec NewCoefImp_trans(const arma::vec& linpredz, const arma::vec& ytrans, const arma::mat& xwxxw){

  arma::vec phiNew(xwxxw.n_rows);

  phiNew = xwxxw * (ytrans - linpredz);
  return( phiNew );
}

double NewLambdaImp(const arma::mat& REs, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const arma::vec& weights, const arma::vec& weightsImp, double tol, int MaxIt, double lambda, const arma::mat& xwxxw){
  arma::uword ind;
  arma::uvec ind2;
  double lambdaNew;
  double lambdaOld = lambda;
  double Hess;
  double GradOld;
  double dir;

  double eps = 1.1 * tol;

  bool klappt;
  List ret;
  double llAlt;
  double ll;
  double llzw;

  arma::uword N = y.n_elem;
  arma::uword m, MI = REs.n_rows;
  double sigmahat;
  arma::vec betahat(X.n_cols);
  arma::vec ytrans(N);
  arma::vec linpredz(N);
  arma::vec res(N);
  bool weightsIn = (weights.n_elem > 1);

  llAlt = loglik_lambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaOld, xwxxw);
  for(int a = 0; a < MaxIt; a++) {

    Hess = HessianLambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaOld, xwxxw);
    GradOld = GradientLambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaOld, xwxxw);

    dir = GradOld/Hess;

    lambdaNew = lambdaOld - dir;
    ll = loglik_lambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaNew, xwxxw);

    eps = abs(dir)/(1+abs(lambdaOld) );

    if(ll > llAlt){
      lambdaNew = lambdaOld - 0.01*dir;

      ll = loglik_lambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaNew, xwxxw);

      if(ll > llAlt){

        lambdaNew = lambdaOld + 0.01*dir;
        ll = loglik_lambdaImp(REs, transform, X, y, Z, weights, weightsImp, lambdaNew, xwxxw);
        if(ll > llAlt) break;
      }
    }
    if(ll > llAlt) break;
    llAlt = ll;
    lambdaOld = lambdaNew;

    if( tol > eps){
      break;
    }

  }

  if(eps > tol) {
    Rcpp::Rcout << "Warning: No convergence for Lambda!" << endl;
    klappt = false;
  }

  ret("phi") = lambdaOld;
  ret("convergence") = klappt;

  return( lambdaOld );
}

List NewCovREImp_trans( const arma::mat& REs, const arma::uvec& qvec, const arma::uvec& ni, const arma::vec& weightsImp){
  arma::uword MI = weightsImp.n_rows;
  arma::uword ns;
  arma::uword q;
  arma::uword b = 0;
  int a = qvec.n_elem;
  arma::mat N;
  arma::mat us;
  bool tries;

  List newCovRE(a);

  for(arma::uword i = 0; i < a; i++) {
    q = qvec.at(i);
    ns = ni.at(i);
    N.set_size(q,q);
    N.zeros();

    us.set_size(q, MI);

    double* mem;
    int usRows = q;

    for(arma::uword j = 0; j < ns; j++) {
      us = REs.cols(b, b+q-1).t();

      for(arma::uword m = 0; m < MI; m++) {
        mem = us.colptr(m);
        arma::vec usCol(mem, usRows, false, true);
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

double evalhv_trans( const arma::vec& REs, const List& invmat, const arma::uvec& qvec, const arma::uvec& ni){
  arma::uword a = qvec.n_elem;
  arma::uword b = 0;
  arma::uword q, ns;
  arma::uword i, j;
  double ret = 0.;

  for(i = 0; i < a; i++){
    q = qvec.at(i);
    ns = ni.at(i);
    arma::mat M = as<arma::mat>(invmat[i]);
    arma::vec us(q);

    for(j = 0; j < ns; j++){
      us = REs.subvec(b, b+q-1);
      ret -= 0.5*as_scalar( us.t()*M*us );
      b += q;
    }
  }
  return(ret);
}

arma::mat ImportanceSampling_trans( int MI, const arma::vec& modus, const arma::mat& Chol, arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale, double lambda, const arma::mat& xwxxw){

  arma::uword nr = Z.n_cols;
  arma::mat rands(nr, MI);
  List control;
  arma::vec re(nr);
  arma::vec us(nr);
  arma::vec w(MI, fill::zeros);
  double sup;
  double sign;
  double addprop;

  arma::uword q;
  arma::uword a = qvec.n_elem;
  List invmat(a);
  double hv;

  double ldet = 0., sgn;

  for(arma::uword i = 0; i < a; i++){
    q = qvec.at(i);
    arma::mat M = as<arma::mat>(CovRE[i]);
    if( q > 1){
      invmat(i) = arma::inv_sympd(M);
    }else{
      invmat(i) = 1/M;
    }

  }

  arma::vec yt = trans(y, transform, lambda);

  for(arma::uword k = 0; k < MI; k = k+2) {
    us.randn();
    addprop = 0.5*as_scalar( us.t() * us ); // Wenn V inverse der Kovarianzmatrix und L*L.t() = V^{-1}, dann ist re.t() * L.t() * V * L * re = re.t() * re
    re = Chol * us;
    arma::vec randsCol = modus + re;
    hv = evalhv_trans( randsCol, invmat, qvec, ni);
    rands.col(k) = randsCol;

    w.at(k) += addprop; // Teile durch proposal dichte
    w.at(k) += hv;
    // w.at(k) -= loglik_lambda(Z * randsCol, transform, X, y, weights, lambda, xwxxw);
    w.at(k) += loglikCond_trans( randsCol, phi, X, yt, Z, weights, scale);

    randsCol = modus - re;
    hv = evalhv_trans( randsCol, invmat, qvec, ni);
    rands.col(k+1) = randsCol;
    w.at(k+1) += addprop;
    w.at(k+1) += hv;
    // w.at(k+1) -= loglik_lambda(Z * randsCol, transform, X, y, weights, lambda, xwxxw);
    w.at(k+1) += loglikCond_trans( randsCol, phi, X, yt, Z, weights, scale);
  }

  rands = rands.t();

  w -= mean(w);

  for( arma::uword i = 0; i < w.n_elem; i++) {
    double value = w.at(i);
    if (value > 35) {
      // Rcpp::Rcout << "A" << endl;
      w.at(i) = 35;
    } else if (value < -25) {
      w.at(i) = -25;
      // Rcpp::Rcout << "B" << endl;
    }
  }

  w = exp(w);
  w /= sum(w);

  rands = join_horiz(rands, w);
  return(rands);

}

List GenModusNeu_trans(const arma::vec& REs, const arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z, const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double tol, int MaxIt, arma::mat Chol, double scale, double lambda, const arma::mat& xwxxw){

  arma::uword nr = REs.n_elem;
  arma::uword ind;
  arma::uvec ind2;
  arma::vec modusNeu(nr);
  arma::vec modusAlt(nr);
  arma::vec gradAlt(nr);
  arma::vec gradNeu(nr);
  arma::vec dir(nr);
  arma::vec res(nr);
  arma::vec alpha = linspace<vec>(0.001,1,25);
  arma::vec ll(25);
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
  arma::vec ytrans = trans(y, transform, lambda);

  List invmat(qvec.n_elem);

  for(arma::uword i = 0; i < qvec.n_elem; i++){
    arma::uword q = qvec.at(i);
    arma::mat M = as<arma::mat>(CovRE[i]);
    if( q > 1){
      invmat[i] = arma::inv_sympd(M);
    }else{
      invmat[i] = 1/M;
    }

  }

  // Special vars
  arma::vec beta = phi;
  double newscale = scale;
  arma::vec re = modusAlt;

  bool weightsIn = weights.n_elem > 1;
  arma::uword n = X.n_rows;
  arma::vec linpredz(n);

  // ende special vars
  // linpredz = Z * re;
  // beta = xwxxw * (ytrans - linpredz);
  // if(weightsIn){
  //   newscale = accu( weights % square(ytrans - X * beta - linpredz))/n;
  // }else{
  //   newscale = accu( square(ytrans - X * beta - linpredz))/n;
  // }

  gradAlt = GradientAll_trans(modusAlt, beta, X, ytrans, Z, CovRE, weights, qvec, ni, newscale, xwxxw);
  dir = -Chol * Chol.t() * gradAlt;
  LLMin = -evalhv_trans(re, invmat, qvec, ni) - loglikCond_trans(re, beta, X, ytrans, Z, weights, scale);
  // LLMin = -evalhv_trans(re, invmat, qvec, ni) + loglik_lambda(linpredz, transform, X, y, weights, lambda, xwxxw);

  double* mem;
  int CholRows = Chol.n_rows;

  for(k = 0; k < MaxIt; k++) {
    for(arma::uword i = 0; i < 25; i++){
      re = modusAlt + dir*alpha.at(i);
      linpredz = Z * re;

      // ll.at(i) = -evalhv_trans(re, invmat, qvec, ni) + loglik_lambda(linpredz, transform, X, y, weights, lambda, xwxxw);
      ll.at(i) = -evalhv_trans(re, invmat, qvec, ni) - loglikCond_trans(re, beta, X, ytrans, Z, weights, scale);
    }

    if(LLMin <= min(ll) || ll.has_nan()){
      break;
    }

    LLMin = min(ll);
    ind2 = sort_index(ll);
    ind = ind2(0);

    dir = dir* alpha.at(ind);
    modusNeu = modusAlt + dir;

    // ############################
    // linpredz = Z * modusNeu;
    // beta = xwxxw * (ytrans - linpredz);
    // if(weightsIn){
    //   newscale = accu( weights % square(ytrans - X * beta - linpredz))/n;
    // }else{
    //   newscale = accu( square(ytrans - X * beta - linpredz))/n;
    // }
    // #########################

    gradNeu = GradientAll_trans(modusNeu, beta, X, ytrans, Z, CovRE, weights, qvec, ni, newscale, xwxxw);
    res = gradNeu - gradAlt;

    if( max(abs(dir)/(1 + abs(modusAlt))) < tol ) {
      conv = true;
      break;
    }
    factor = 1/accu( res % dir);

    s = Chol.t() * res;
    ind2 = find(s, 1, "last");
    ind = ind2(0L);

    h = s.at(ind) * Chol.col(ind);
    phifac = s.at(ind)*s.at(ind);

    Chol.col(0L) = dir * sqrt(factor);
    sigmaNeu = min( sigmaAlt, norm(Chol.col(0L)) );

    for(arma::uword i = ind; i > 0; i--) {
      mem = Chol.colptr(i-1);
      arma::vec CholCol(mem, CholRows, false, false);
      Chol.col(i) = sqrt(phifac/(s.at(i-1)*s.at(i-1) + phifac)) * ( s.at(i-1)/phifac * h - CholCol);
      // Chol.col(i) = Chol.col(i) - accu(res % CholCol) * factor * dir;
      Chol.col(i) *= sigmaNeu / sqrt(accu(square(Chol.col(i))));
      phifac += s(i-1)*s.at(i-1);
      h += s.at(i-1)*CholCol;
    }

    sigmaAlt = sigmaNeu;
    // Hier die Cholesky der INVERSEN der Hesse-Matrix
    modusAlt = modusNeu;
    gradAlt = gradNeu;
    dir = -Chol*Chol.t()*gradAlt;
  }

  if(!conv) {
    Rcpp::Rcout << "No convergence!" << endl;
  }
  ret("modus") = modusAlt;
  ret("invChol") = Chol;
  ret("iter") = k;

  return(ret);
}

// ########################################################################################################################################################
// ########################################################################################################################################################
// ########################################################################################################################################################
// Likelihood (not expected LL)
double likImp(int MI, const arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni, double scale, double lambda){
  arma::uword a = qvec.n_elem;
  arma::uword nr = Z.n_cols;
  arma::uword b = 0;
  arma::uword i, j, m, q, ns;
  arma::vec us;
  arma::mat REs(nr,MI);
  arma::mat M;
  double likeli = 0.;
  double adj = 0.;

  arma::vec yt = trans(y, transform, lambda);
  double scaling = loglikCond_trans(zeros(nr), phi, X, yt, Z, weights, scale);
  double recomp = 0.;

  for( i = 0; i < a; i++){
    q = qvec.at(i);
    ns = ni.at(i);
    us.set_size(q);
    M = as<arma::mat>(CovRE[i]);
    M = chol(M);

    recomp -= 2 * ns * log(arma::accu(M.diag()));

    for(j = 0; j < ns; j++){
      for(m = 0; m < MI; m=m+2){
        us.randn();
        likeli += exp(-accu( square(us)) );
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
    double ll = loglikCond_trans(us, phi, X, yt, Z, weights, scale);
    likeli += exp(ll - scaling);
  }
  likeli /= MI;

  likeli = log(likeli) + recomp + scaling;

  if(transform == "box-cox"){
    if(weights.n_elem == 1) adj = (lambda - 1) * arma::accu( log(y));
    if(weights.n_elem > 1) adj = (lambda - 1) * arma::accu( weights % log(y));
  }else if(transform == "dual"){
    arma::vec ypow = pow(y, lambda);
    if(weights.n_elem == 1) adj = (lambda - 1) * arma::accu( log(ypow + 1/ypow) - log(y));
    if(weights.n_elem > 1) adj = (lambda - 1) * arma::accu( weights % (log(ypow + 1/ypow) - log(y)) );
  }

  likeli += adj;

  return(likeli);
}

// #####################################################################################################
arma::mat constrBigInvCov(const List& CovRE, const arma::uvec& qvec, const arma::uvec& ni){
  arma::uword nr = accu( qvec % ni);
  arma::uword q = qvec.n_elem, i, ns, qs, b = 0, j;
  arma::mat M;
  arma::mat ret(nr, nr, fill::zeros);
  for(i = 0; i < q; i++){
    qs = qvec(i);
    ns = ni(i);
    M.set_size(qs, qs);
    M = as<arma::mat>(CovRE[i]);
    if( qs > 1) M = arma::inv_sympd(M); else M = 1/M;
    // if(qs > 1) M = chol(M); else M = sqrt(M);
    for(j = 0; j < ns; j++){
      ret.submat(b, b, b+qs-1, b+qs-1) = M;
      b += qs;
    }
  }
  return(ret);
}

arma::mat constrBigCov(const List& CovRE, const arma::uvec& qvec, const arma::uvec& ni){
  arma::uword nr = accu( qvec % ni);
  arma::uword q = qvec.n_elem, i, ns, qs, b = 0, j;
  arma::mat M;
  arma::mat ret(nr, nr, fill::zeros);
  for(i = 0; i < q; i++){
    qs = qvec(i);
    ns = ni(i);
    M.set_size(qs, qs);
    M = as<arma::mat>(CovRE[i]);
    for(j = 0; j < ns; j++){
      ret.submat(b, b, b+qs-1, b+qs-1) = M;
      b += qs;
    }
  }
  return(ret);
}

arma::mat constrBigChol(const List& CovRE, const arma::uvec& qvec, const arma::uvec& ni){
  arma::uword nr = accu( qvec % ni);
  arma::uword q = qvec.n_elem, i, ns, qs, b = 0, j;
  arma::mat M;
  arma::mat ret(nr, nr, fill::zeros);
  for(i = 0; i < q; i++){
    qs = qvec(i);
    ns = ni(i);
    M.set_size(qs, qs);
    M = as<arma::mat>(CovRE[i]);
    M = chol(M, "lower");
    for(j = 0; j < ns; j++){
      ret.submat(b, b, b+qs-1, b+qs-1) = M;
      b += qs;
    }
  }
  return(ret);
}

// ###############################################################################
// [[Rcpp::export(.wlmm_transform_cpp)]]
List wlmm_transform_cpp( const arma::vec& phi, const String& transform, arma::mat& X, const arma::vec& y, const arma::sp_mat& Z,  const List& CovRE, const arma::vec& weights, const arma::uvec& qvec, const arma::uvec& ni,  double scale, double lambda, int MI, int iter1, int iter2, double tol1, double tol2, bool trace, int nDecrease){

  List ret;
  List outModus;
  List covmat;
  List covmatAlt;
  List covmatAllTimeBest;

  arma::uword nq = Z.n_cols;
  arma::uword nqm1 = nq-1L;
  arma::uword ngrps = qvec.n_elem;
  arma::uword counter1 = 0;
  arma::uword counter3 = 0;
  arma::uword p = phi.n_elem;
  arma::uword n = y.n_elem;
  bool weightsIn = (weights.n_elem > 1);

  double eps1 = 1.1*tol1;
  double scaleAlt;
  double scaleNeu;
  double scaleAllTimeBest;
  double lambdaAlt;
  double lambdaNeu;
  double lambdaAllTimeBest;
  double ll;
  double llAlt;
  double llAllTimeBest;

  arma::vec phiNeu(p);
  arma::vec phiAlt(p);
  arma::vec hv(p+1);
  arma::vec phiAllTimeBest(p);
  arma::vec modus(nq, fill::zeros);
  arma::vec ytrans(n);
  arma::vec ytransAlt(n);
  arma::vec res(n);
  arma::vec pred(n);

  arma::mat xwxxw(X.n_cols, X.n_rows);
  // arma::mat cholinvZWZ(nq, nq);
  if(weightsIn){
    xwxxw = arma::inv_sympd( X.t() * arma::diagmat(weights) * X ) * X.t() * arma::diagmat(weights);
    // cholinvZWZ = Z.t() * arma::diagmat(weights) * Z;
  }else{
    xwxxw = arma::inv_sympd( X.t() * X ) * X.t();
    // cholinvZWZ = Z.t() * Z;
  }

  // cholinvZWZ = chol( arma::inv_sympd(cholinvZWZ) );

  arma::mat Us(MI, nq + 1L);
  arma::mat Chol = eye<arma::mat>(nq,nq);

  Chol = constrBigChol(CovRE, qvec, ni);
  covmatAlt = CovRE;
  covmatAllTimeBest = CovRE;

  phiAlt = phi;
  phiAllTimeBest = phi;
  scaleAlt = scale;
  scaleAllTimeBest = scale;
  lambdaAlt = lambda;
  lambdaAllTimeBest = lambda;

  int mi = 0.1 * MI;
  if( (mi/2)*2 != mi ){
    mi++;
  }
  bool updateLambda = (transform != "none");

  ytransAlt = trans(y, transform, lambdaAlt);

  Chol = eye(nq,nq) /norm(GradientAll_trans(modus, phiAlt, X, ytransAlt, Z, covmatAlt, weights, qvec, ni, scaleAlt, xwxxw));

  outModus = GenModusNeu_trans( modus, phiAlt, transform, X, y, Z, covmatAlt, weights, qvec, ni, tol2, iter2, Chol, scaleAlt, lambdaAlt, xwxxw);
  Chol = as<arma::mat>(outModus[1]);
  modus = as<arma::vec>(outModus[0]);

  Us = ImportanceSampling_trans( MI, modus, Chol, phiAlt, transform, X, y, Z, covmatAlt, weights, qvec, ni, scaleAlt, lambdaAlt, xwxxw);
  llAlt = logliktTot_trans(modus, phiAlt, transform, X, y, Z, covmatAlt, weights, qvec, ni, scaleAlt, lambdaAlt, xwxxw);
  llAllTimeBest = llAlt;

  while(counter1 < iter1 & eps1 > tol1){
    int counter2 = 0;
    for(int counter2 = 0; counter2 < 10; counter2++ ){

      Us = ImportanceSampling_trans( MI + counter2 * mi, modus, Chol, phiAlt, transform, X, y, Z, covmatAlt, weights, qvec, ni, scaleAlt, lambdaAlt, xwxxw);
      pred.zeros();

      for(arma::uword m = 0; m < MI + counter2 * mi; m++){
        pred += Z * Us.row(m).head(nq).t() * Us(m, nq);
      }

      if(updateLambda) lambdaNeu = NewLambdaImp( Us.cols(0,nqm1), transform, X, y, Z, weights, Us.col(nq), tol2, iter2, lambdaAlt, xwxxw);
      ytrans = trans(y, transform, lambdaNeu);

      covmat = NewCovREImp_trans( Us.cols(0, nqm1), qvec, ni, Us.col(nq));
      phiNeu = NewCoefImp_trans( pred, ytrans, xwxxw);
      if(phiNeu.has_nan()){
        break;
      }
      res = ytrans - X * phiNeu;
      arma::vec linpredz = Z * sum(Us.cols(0,nqm1).each_col() % Us.col(nq), 0).t();
      scaleNeu = 0;
      if(weightsIn){
        for(arma::uword m = 0; m < Us.n_rows; m++){
          arma::vec us = Us.row(m).t();
          scaleNeu +=  us(nq) * arma::accu(weights % arma::square(res - Z * us.head(nq)) );
        }
        scaleNeu /= accu(weights);

      }else{
        for(arma::uword m = 0; m < Us.n_rows; m++){
          arma::vec us = Us.row(m).t();
          scaleNeu +=  us(nq) * arma::accu( arma::square(res - Z * us.head(nq)) );
        }
        scaleNeu /= n;
      }

      Us = ImportanceSampling_trans( MI + counter2 * mi, modus, Chol, phiAlt, transform, X, y, Z, covmatAlt, weights, qvec, ni, scaleAlt, lambdaAlt, xwxxw);
      ll = logliktTotImp_trans(Us.cols(0,nqm1), phiNeu, transform, X, y, Z, covmat, weights, qvec, ni, Us.col(nq), scaleNeu, lambdaNeu, xwxxw);

      if( ll < llAlt | counter1 == 0){
        break;
      }
    }
    if(phiNeu.has_nan()){
      break;
    }
    ytransAlt = ytrans;

    // modus = zeros(nq);
    modus = sum(Us.cols(0,nqm1).each_col() % Us.col(nq)).t();
    // Chol = constrBigCov(covmat, qvec, ni);
    Chol = eye(nq,nq) /norm(GradientAll_trans(modus, phiNeu, X, ytrans, Z, covmat, weights, qvec, ni, scale, xwxxw));
    outModus = GenModusNeu_trans( modus, phiNeu, transform, X, y, Z, covmat, weights, qvec, ni, tol2, iter2, Chol, scaleNeu, lambdaNeu, xwxxw);
    Chol = as<arma::mat>(outModus[1]);
    modus = as<vec>(outModus[0]);

    MI += mi;

    eps1 = max( arma::vectorise( abs(as<arma::mat>(covmat[0]) - as<arma::mat>(covmatAlt[0])) / ( abs(as<arma::mat>(covmatAlt[0]))  + 1) ) );

    for(int i = 0; i< ngrps;i++){
      double eps = max( arma::vectorise( abs(as<arma::mat>(covmat[i]) - as<arma::mat>(covmatAlt[i])) / ( abs( as<arma::mat>(covmatAlt[i]) )  + 1)  ) ) ;
      if(eps > eps1){
        eps1 = eps;
      }
    }

    eps1 = std::max( max( abs(phiNeu - phiAlt) / ( abs(phiAlt) + 1) ), eps1);
    eps1 = std::max( eps1, abs(lambdaNeu-lambdaAlt)/(1+abs(lambdaAlt)) );
    eps1 = std::max( eps1, abs(scaleNeu - scaleAlt)/(1 + abs(scaleAlt)) );

    if(ll < llAlt | counter1 == 0){
      llAllTimeBest = ll;
      phiAllTimeBest= phiNeu;
      covmatAllTimeBest = covmat;
      scaleAllTimeBest = scaleNeu;
      if(updateLambda) lambdaAllTimeBest = lambdaNeu;
    }else{
      counter3++;
    }

    if(counter3 == nDecrease){
      Rcpp::Rcout << "Keine Verbesserung der Likelihood mehr!" << endl;
      break;
    }

    if( (mi/2)*2 != mi ){
      mi++;
    }
    Us.set_size(MI,nq+1);

    counter1++;
    if(trace){
      Rcpp::Rcout << "Iteration " << counter1 << " von " << iter1<< endl;
      Rcpp::Rcout << "Veraenderung zum letzten Schritt: " << eps1 << endl;
      Rcpp::Rcout << phiNeu << endl;
      Rcpp::Rcout << "Lambda: " << lambdaNeu << endl;
      for(int l = 0; l < qvec.size(); l++){
        Rcpp::Rcout << as<arma::mat>(covmat[l]) << endl;
      }
      Rcpp::Rcout << "Residual Variance: " << scaleNeu << endl;
      Rcpp::Rcout << "Increased Log-Likelihood: " << (ll < llAlt) << endl;
    }

    llAlt = ll;
    phiAlt = phiNeu;
    scaleAlt = scaleNeu;
    covmatAlt = covmat;

    if(updateLambda) lambdaAlt = lambdaNeu;
    ytrans = trans(y, transform, lambdaAlt);

  }

  ytrans = trans(y, transform, lambdaAllTimeBest);

  // modus = zeros(nq);
  // Chol = constrBigChol(covmatAllTimeBest, qvec, ni);
  modus = sum( Us.cols(0,nqm1).each_col() % Us.col(nq) ).t();
  Chol = 1/ sqrt( norm(GradientAll_trans(modus, phiAllTimeBest, X, ytrans, Z, covmatAllTimeBest, weights, qvec, ni, scaleAllTimeBest, xwxxw)) ) * eye(nq,nq);

  outModus = GenModusNeu_trans( modus, phiAllTimeBest, transform, X, y, Z, covmatAllTimeBest, weights, qvec, ni, tol2, iter2, Chol , scaleAllTimeBest, lambdaAllTimeBest, xwxxw);
  Chol = as<arma::mat>(outModus[1]);
  modus = as<vec>(outModus[0]);

  Us = ImportanceSampling_trans( MI + mi, modus, Chol , phiAllTimeBest, transform, X, y, Z, covmatAllTimeBest, weights, qvec, ni, scaleAllTimeBest, lambdaAllTimeBest, xwxxw);

  ll = logliktTotImp_trans(Us.cols(0,nqm1), phiAllTimeBest, transform, X, y, Z, covmatAllTimeBest, weights, qvec, ni, Us.col(nq), scaleAllTimeBest, lambdaAllTimeBest, xwxxw);
  arma::vec one = ones(1,1);
  double ll_mode = logliktTot_trans(modus, phiAllTimeBest, transform, X, y, Z, covmatAllTimeBest, weights, qvec, ni, scaleAllTimeBest, lambdaAllTimeBest, xwxxw);

  // Berechnung der Kovarianzmatrix der Fixed Effects:
  arma::mat gr2_beta(phi.n_elem, phi.n_elem, fill::zeros);
  arma::mat Hess_beta = HessianBetaImp_trans( X, weights, scaleAllTimeBest);
  arma::vec gr_beta = GradientBetaImp_trans(Us.cols(0,nqm1), phiAllTimeBest, X, ytrans, Z, weights, Us.col(nq), scaleAllTimeBest);

  for(arma::uword b = 0; b < Us.n_rows; b++){
    rowvec re = Us.row(b);
    arma::vec gr_b = GradientBetaImp_trans( re.head(nq), phiAllTimeBest, X, ytrans, Z, weights, one, scaleAllTimeBest);
    gr2_beta += gr_b * gr_b.t() * as_scalar(re.tail(1)) ;
  }
  arma::mat fisher_beta = Hess_beta - gr2_beta + gr_beta * gr_beta.t();
  fisher_beta = inv(fisher_beta);

  ret("coef") = phiAllTimeBest;
  ret("vcov_coef") = fisher_beta;
  ret("VarCov") = covmatAllTimeBest;
  ret("lambda") = lambdaAllTimeBest;
  ret("sigma") = sqrt( scaleAllTimeBest );
  ret("transform") = transform;
  ret("RE_mat") = Us;
  ret("RE_mode") = modus;
  ret("RE_invHessian") = Chol * Chol.t();
  ret("LLmod") = -ll_mode;
  ret("LLexp") = -ll;
  ret("LL") = likImp(MI, phiAllTimeBest, transform, X, y, Z,  covmatAllTimeBest, weights, qvec, ni, scaleAllTimeBest, lambdaAllTimeBest);
  ret("convergence") = eps1 < tol1;

  return(ret);
}












