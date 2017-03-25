#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List compressCpp(NumericMatrix X, int q, NumericVector Y, int s) {
  int p = X.ncol();
  int n = X.nrow();

  NumericMatrix QX(q, p);
  NumericVector QY(q);
  double u1, u2;
  double rescale = sqrt(q/s);
  double prob = 1.0 - 1.0/s;

  RNGScope scope;
  for(int qiter = 0; qiter < q; qiter++){
    for(int niter = 0; niter < n; niter++){
      u1 = runif(1)[0];
      if(u1 < prob) continue;
      u2 = runif(1)[0];
      if(u2 < 0.5){
        QY[qiter] += Y[niter];
        for(int jiter = 0; jiter < p; jiter++){
          QX(qiter, jiter) += X(niter, jiter);
        }
      }else{
        QY[qiter] -= Y[niter];
        for(int jiter = 0; jiter < p; jiter++){
          QX(qiter, jiter) -= X(niter, jiter);
        }
      }
    }
  }

  return List::create(Named("QX") = QX / rescale,
                      Named("QY") = QY / rescale);
}


