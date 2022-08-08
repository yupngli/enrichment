#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double predprobPosterior(int x, int n, int nmax, double alpha, double beta, double p0, double theta_t){
  
  double prob = 0.0;
  double eps = std::numeric_limits<double>::epsilon();
  double pxy;
  
  for (int resp = 0; resp < nmax - n + 1; resp++) {
    pxy = (1.0 - R::pbeta(p0, alpha + resp + x, beta + nmax - resp - x, 1, 0));
    if (pxy > theta_t || std::abs(pxy - theta_t) < eps) {
      prob += exp(
        R::lchoose(nmax - n, resp) +
          R::lbeta(alpha + resp + x, beta + nmax - resp - x) -
          R::lbeta(alpha + x, beta + n - x)
      );
    }
  }
  return prob;
}

