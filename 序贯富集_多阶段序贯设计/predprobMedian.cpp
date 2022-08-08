#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double predprobMedian(int x, int n, int nmax, double alpha, double beta, double theta_t) {

  double prob = 0.0;
  double eps = std::numeric_limits<double>::epsilon();
  double med;

  for (int resp = 0; resp < nmax - n + 1; resp++) {
    // --------------------------------------------------------
    // --- use median instead of posterior in pfizer's settings
    // --------------------------------------------------------
    // generate beta random numbers
    NumericVector rand = Rcpp::rbeta( 100000, alpha + x + resp, beta + nmax - x - resp);
    med = median(rand);
    if (med > theta_t || std::abs(med - theta_t) < eps) {
      prob += exp(
        R::lchoose(nmax - n, resp) +
          R::lbeta(alpha + x + resp, beta + nmax - x - resp) -
          R::lbeta(alpha + x, beta + n - x)
      );
    }
  }

  return prob;

}


