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

void r_sbern(arma::rowvec &Q, double s, double rs){
  double p = 1.0 / s;
  double q = 1.0 - p;
  
  arma::rowvec::iterator it = Q.begin();
  arma::rowvec::iterator it_end = Q.end();
  
  for(; it < it_end; it++){
    double u = rng_unif();
    if(u < p) (*it) += rs;
    if(u > q) (*it) -= rs;
  }
}

// [[Rcpp::export]]
List compressCpp(arma::mat const & X, int q, arma::colvec const & y, double s){
  int p = X.n_cols;
  int n = X.n_rows;
  
  arma::mat QX(q, p, arma::fill::zeros);
  arma::vec QY(q, arma::fill::zeros);
  
  double rescale = sqrt(s/(1.0*q));
  
  for(int i=0; i<q; i++){
    arma::rowvec Q(n, arma::fill::zeros);
    r_sbern(Q, s, rescale);
    QX.row(i) += Q * X;
    QY(i) += arma::as_scalar(Q * y);
  }
  return List::create(Named("QX") = QX, Named("QY") = QY);
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
