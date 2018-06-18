
#include "general_util/sigfig.h"
#include "calcem/calcem_util/sig.h"
#include <cmath>

auto do150( int& i, std::vector<double>& x, std::vector<double>& y, 
  double& xm, double& ym, double& yt, double& test, const double& tolmin, 
  const double& e, const double& u, const double& tev, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<std::vector<double>>& sab, const double& tevz, const int& lasym,
  const double& az, const double& az2, const double& teff2, const int lat, const double& cliq,
  const double& sb, const double& sb2, const double& teff, const double& tol, int iinc){
  bool continue_to_150 = false;
  //std::cout << std::setprecision(20) << 150 << "     " << y[0] << std::endl;
  if (i <= 3 or 0.5*(y[i-2]+y[i-1])*(x[i-2]-x[i-1]) >= tolmin) {
    xm = 0.5*(x[i-2]+x[i-1]);
    xm = sigfig(xm,8,0);
    if (xm > x[i-1] and xm < x[i-2]){
      ym=0.5*(y[i-2]+y[i-1]);
      yt = sig( e, xm, u, tev, alpha, beta, sab, az, tevz, lasym, 
                az2, teff2, lat, cliq, sb, sb2, teff, iinc );
      test = tol*std::abs(yt);
      
      if (std::abs(yt-ym) > test) {
        // point fails
        i=i+1;
        x[i-1]=x[i-2];
        y[i-1]=y[i-2];
        x[i-2]=xm;
        y[i-2]=yt;
        // go to 150
        // continue;
        continue_to_150 = true;
      }
    }
  }
  return continue_to_150;
}
 

