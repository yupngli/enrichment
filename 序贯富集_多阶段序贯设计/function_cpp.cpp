#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector responsegate(int maxresponse,int samplesizestart,int maxsamplesize,double a,double b,double nv,double positivestatsig,double statsig,double postmed){
  // Creating object
  int row,totalrow,f1,f2,f3;
  double med,sig;
  totalrow = (maxsamplesize-samplesizestart+1)*(maxresponse);
  NumericMatrix df(totalrow,8);
  row = 0;
  // using loop to construct response matrix
  for(int r = 1; r <= maxresponse; ++r) {
    for(int n = samplesizestart; n <= maxsamplesize; ++n) {
      // generate beta random numbers
      NumericVector rand = Rcpp::rbeta( 100000, a+r, b+n-r);
      // compare to nv
      LogicalVector logicnv = rand >= nv;
      sig= mean(logicnv);
      med = median(rand);
      // output
      df(row,0) = r;
      df(row,1) = n;
      df(row,2)=sig;
      df(row,3)=med;
      // decision flag
      f1=0;
      f2=0;
      f3=0;
      if (sig<positivestatsig && med>=postmed){
        df(row,4)=1;
        f1=1;
        };
      if (sig>=positivestatsig && med>=postmed && sig<statsig){
        df(row,5)=1;
        f2=1;
        };
      if (sig>=statsig && med>=postmed){
        df(row,6)=1;
        f3=1;
        };
      // 2-full,3-positive only,4-n/a
      if (f3==1){
        df(row,7) = 2;
      } else if(f2==1){
        df(row,7) = 3;
      } else if (f1==1){
        df(row,7) = 4;
      }
      row++;
    }
  }
  return df;
}
