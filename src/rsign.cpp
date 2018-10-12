#include <Rcpp.h>
using namespace Rcpp;


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
