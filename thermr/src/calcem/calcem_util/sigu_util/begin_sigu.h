#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"

auto do_113_116( int& jbeta, const int& lat, std::vector<double>& x, 
  std::vector<double>& y, const double& e, const double& tev, const double& tevz,
  const double& root1, const double& u,
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<std::vector<double>>& sab, const double& az, 
  const int lasym, const double& az2, const double& teff, const double& teff2, 
  const double& cliq, const double& sb, const double& sb2, const int& iinc){

    std::cout << 113 << std::endl;
  do {
    // 113 continue
    if (jbeta == 0) jbeta=1;
    if (lat == 1) { x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tevz; }
    else {          x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tev;  }
    x[0] = sigfig(x[0],8,0);
    if ( jbeta <= 0 and x[0] == e ){ x[0] = sigfig(e,8,-1); }
    ++jbeta;

  } while(x[0] <= x[1]);

  jbeta = jbeta-1;
   
  std::cout << 116 << std::endl;
  if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
    x[0]=root1*root1;
  }

  //x[0] = sigfig(x[0],8,0);
  x[0] = sigfig(x[0],8,1);
     //std::cout << std::setprecision(25) << x[0] << "    " << x[1] << "     " << x[2] << std::endl;
  //std::cout << sb2 << "        " << iinc << std::endl;
  std::cout << std::setprecision(25) << "x[0]    " << x[0] << std::endl;
  y[0] = sig( e, x[0], u, tev, alpha, beta, sab, az, tevz, lasym, az2,
      teff2, lat, cliq, sb, sb2, teff, iinc );
} 



