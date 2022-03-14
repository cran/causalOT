#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double otl2_obj(NumericVector & vars, IntegerVector & f_idx,
                IntegerVector & g_idx, NumericVector & a,
                NumericVector & b, 
                NumericMatrix & cost, double lambda,
                int N, int M) {
  
  NumericMatrix diff(N,M);
  diff.fill(0.0);
  for (int i = 0; i < N; i ++) {
    for ( int j = 0; j < M ; j++) {
      double calc = vars(f_idx(i)) + vars(g_idx(j)) - cost(i,j);
      if (calc < 0.0) {
        calc = 0.0;
      } else {
        calc = calc * calc;
      }
      diff(i,j) = calc;
    }
  }
  
  double obj = 0.0;
  
  for(int i = 0; i < N ; i++) {
    obj += vars(f_idx(i)) * a(i);
  }
  
  for(int j = 0; j < M ; j++) {
    obj += vars(g_idx(j)) * b(j);
  }
  
  obj -= 0.5 / lambda * sum(diff);
  return(obj);
  
}
