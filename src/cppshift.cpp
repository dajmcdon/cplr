#include <RcppArmadillo.h>
#include <cstdlib> 
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

// [[Rcpp::export]]
NumericVector r_sign(const int& n) {
  NumericVector x(n);
  
  for (int i = 0; i < n; i++){
    double u = rng_unif();
    x[i] = (u > 0.5) ? 1.0 : -1.0;
  }
  return x;
}


//' @export
// [[Rcpp::export]]
arma::mat getDesign(arma::cube const& arr, arma::umat mask, int rad){
  
  arma::cube work(arr);
  mask--; // R indexing issue.
  
  // padding
  work.insert_rows(arr.n_rows-1, rad);
  work.insert_rows(0,rad);
  work.insert_cols(arr.n_cols-1,rad);
  work.insert_cols(0,rad);
  work.insert_slices(arr.n_slices-1,rad);
  work.insert_slices(0,rad);
  // Rcout << "Padded" << std::endl;
  // Rcout << "Work nrow" << work.n_rows << std::endl;
  // Rcout << "Work ncol" << work.n_cols << std::endl;
  // Rcout << "Work nsli" << work.n_slices << std::endl;
  
  // increment mask indices
  mask += rad;
  arma::umat mu(mask);
  arma::umat md(mask);
  mu += rad;
  md -= rad;
  // Rcout << "Re masked" << std::endl;
  
  int s = 2*rad + 1;
  int sz = s*s*s;
  // Rcout << "Out nrow" << mask.n_rows << std::endl;
  // Rcout << "Out ncol" << sz << std::endl;
  arma::mat out(mask.n_rows, sz);
  out.zeros();
  arma::cube wc(s,s,s);
  for(int iter=0; iter<mask.n_rows; iter++){
    wc = work.subcube(md(iter,0),md(iter,1),md(iter,2), 
                      mu(iter,0),mu(iter,1),mu(iter,2));
    out.row(iter) += arma::vectorise(wc).t();
  }
  // kill self
  int self = floor(sz/2);
  out.shed_col(self);
  return out;
}
