#include "general_util/sigfig.h"


auto do313( int& jbeta, const int& lat, double& ep, double& enow, 
  const std::vector<double>& beta, const double& tev, const double& tevz,
  const std::vector<double>& x, int& iskip){

  double tev_val = lat == 1 ? tevz : tev;
  do {
    //std::cout << 313 << std::endl;
    if (jbeta == 0) jbeta = 1;
    if (jbeta < 0) {
      ep = enow - beta[-jbeta-1] * tev_val;
      ep = (ep == enow) ? sigfig(enow,8,-1) : sigfig(ep,8,0);
    } 
    else {
      ep = enow + beta[jbeta-1] * tev_val;
      if   (ep == enow)   { iskip = 1; }
      ep = (ep == enow) ? sigfig(enow,8,+1) : sigfig(ep,8,0);

    } // endif
    ++jbeta;
  } while (ep <= x[1]);

  jbeta -=1;

}


