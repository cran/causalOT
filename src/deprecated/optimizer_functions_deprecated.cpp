// [[Rcpp::export]]
double cotDual_obj_(const SEXP & vars_, 
                    const SEXP & b_,
                    const SEXP & delta_,
                    const SEXP & BalTarg_,
                    const SEXP & BalMat_,
                    const SEXP& cost_, 
                    double lambda,
                    const std::string penalty) {
  // const SEXP & delta_pos_,
  // const SEXP & delta_neg_,
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  // const vecMap delta_pos(Rcpp::as<vecMap>(delta_pos_));
  // const vecMap delta_neg(Rcpp::as<vecMap>(delta_neg_));
  const vecMap delta(Rcpp::as<vecMap>(delta_));
  const vecMap BalTarg(Rcpp::as<vecMap>(BalTarg_));
  const matMap BalMat(Rcpp::as<matMap>(BalMat_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  //get omega fun
  Rcpp::XPtr<omegaPtr> xpfun = getPtrOmega(penalty);
  Rcpp::XPtr<omegaPtr> xpfun_d = getPtrOmega_d(penalty);
  omegaPtr omega_fun = *xpfun;
  omegaPtr d_omega_fun = *xpfun_d;
  
  // dimensions
  int N = cost.rows();
  int M = cost.cols();
  int K = BalMat.cols();
  
  // assign variables
  vector g = vars.block(0,0,N,1);
  vector xi_pos = vars.block(N,0,K,1);
  vector xi_neg = vars.tail(K);
  
  // xi variable
  vector xi = xi_pos - xi_neg;
  
  // setup output vector
  double output = 0.0;
  
  // Long portion doing xi lasso checks
  // set xi variable to 0 if necessary
  vector xi_grad_check = BalTarg;
  vector g_vec(N);
  
  for(int j = 0; j < N; j ++) {
    g_vec.fill(g(j));
    vector eta_vec = g_vec - BalMat * xi - cost.col(j);
    xi_grad_check += BalMat.transpose() * d_omega_fun(eta_vec, lambda);
  }
  
  vector xi_nzero = (xi_grad_check.cwiseAbs().array() <= delta.array()).cast<double>();
  xi.array() *= xi_nzero.array();
  
  // calculate omega gradient
  for(int j = 0; j < N; j ++) {
    g_vec.fill(g(j));
    vector eta_vec = g_vec - BalMat * xi - cost.col(j);
    vector out_omega = omega_fun(eta_vec, lambda);
    output += out_omega.sum();
  }
  
  // add in penalty parts
  vector sign = xi.cwiseSign();
  
  output += g.dot(b);
  output += -delta.dot(sign) - BalTarg.dot(xi);
  
  return( -output ); // negative b/c lbfgs minimizes by default and we need to maximize
}

// [[Rcpp::export]]
Rcpp::NumericVector cotDual_grad_(const SEXP & vars_, 
                                  const SEXP & b_,
                                  const SEXP & delta_,
                                  const SEXP & BalTarg_,
                                  const SEXP & BalMat_,
                                  const SEXP& cost_, 
                                  double lambda,
                                  const std::string penalty) {
  // const SEXP & delta_pos_,
  // const SEXP & delta_neg_,
  
  // map R types to Eigen types
  const vecMap vars(Rcpp::as<vecMap>(vars_));
  const vecMap b(Rcpp::as<vecMap>(b_));
  // const vecMap delta_pos(Rcpp::as<vecMap>(delta_pos_));
  // const vecMap delta_neg(Rcpp::as<vecMap>(delta_neg_));
  const vecMap delta(Rcpp::as<vecMap>(delta_));
  const vecMap BalTarg(Rcpp::as<vecMap>(BalTarg_));
  const matMap BalMat(Rcpp::as<matMap>(BalMat_));
  const matMap cost(Rcpp::as<matMap>(cost_));
  
  //get omega fun
  Rcpp::XPtr<omegaPtr> xpfun = getPtrOmega(penalty);
  Rcpp::XPtr<omegaPtr> xpfun_d = getPtrOmega_d(penalty);
  omegaPtr omega_fun = *xpfun;
  omegaPtr d_omega_fun = *xpfun_d;
  
  // dimensions
  int N = cost.rows();
  int M = cost.cols();
  int K = BalMat.cols();
  
  // assign variables
  vector g = vars.block(0,0,N,1);
  vector xi_pos = vars.block(N,0,K,1);
  vector xi_neg = vars.tail(K);
  
  // xi variable
  vector xi = xi_pos - xi_neg;
  vector sign = xi.cwiseSign();
  
  // setup output vector
  vector grad(vars.rows());
  vector g_grad = b;
  vector xi_grad = vector::Zero(K);
  
  // Long portion doing xi lasso checks
  // set xi variable to 0 if necessary
  vector xi_grad_check = BalTarg;
  vector g_vec(N);
  
  for(int j = 0; j < N; j ++) {
    g_vec.fill(g(j));
    vector eta_vec = g_vec - BalMat * xi - cost.col(j);
    xi_grad_check += BalMat.transpose() * d_omega_fun(eta_vec, lambda);
  }
  
  vector xi_nzero = (xi_grad_check.cwiseAbs().array() <= delta.array()).cast<double>();
  xi.array() *= xi_nzero.array();
  
  // calculate omega gradient
  for(int j = 0; j < N; j ++) {
    g_vec.fill(g(j));
    vector eta_vec = g_vec - BalMat * xi - cost.col(j);
    vector out_omega = d_omega_fun(eta_vec, lambda);
    g_grad(j) += out_omega.sum();
    xi_grad += BalMat.transpose() * out_omega;
  }
  
  // separate xi pos and neg gradients (also add in penalty terms as appropriate)
  vector xi_grad_pos = xi_grad; // positive portion of xi_grad
  vector xi_grad_neg = -xi_grad; // negative portion of xi_grad
  
  for(int k = 0; k < K; k ++) {
    if (xi(k) < 0.) {
      xi_grad_neg(k) += -BalTarg(k) - delta(k); // only increment negative protion because xi_pos should be zero...
    } else if (xi(k) > 0.) {
      xi_grad_pos(k) += BalTarg(k) + delta(k);
    }
  }
  
  // concatenate gradients
  grad.block(0,0,N,1) = g_grad;
  grad.block(N,0,K,1) = xi_grad_pos;
  grad.block(N+K,0,K,1) = xi_grad_neg;
  
  // convert to R vector
  Rcpp::NumericVector output = Rcpp::wrap(-grad); 
  // negative b/c lbfgs minimizes by default and we need to maximize
  
  return(output);
}
