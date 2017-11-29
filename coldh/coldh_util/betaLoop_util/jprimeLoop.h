#include <iostream>
#include <vector>
#include "../../../discre/discre_util/sint.h"
#include "jprimeLoop_util/sumh.h"

auto jPrime( double& total, int j, const double& be, const double& x, 
  const double& sw, const double& pj, int jj, const std::vector<double>& bex, 
  const std::vector<double>& rdbex, const std::vector<double>& sex, 
  const std::vector<double>& betan, const double& al, const double& wt, 
  const double& tbart, const double& y, int nbx, bool odd ){

  //--sum over the odd or even values of j-prime
  
  double betap, add, snl = 0, pi = 3.14159265358979, tmp, bn;
  int jp, start, end;

  // Select start and end values for loop. If we want odd values, 1-->9,
  // if even 2-->8.
  if ( odd == true ){ start = 1; end = 10; }
  else              { start = 0; end = 9;  }

  for ( auto jp = start; jp < end; jp = jp + 2 ){

    betap = ( -j * (j+1) + jp * (jp+1) ) * x * 0.5;
    bn = be + betap;

    tmp = (2*jp+1) * pj * sw * 4 * sumh(j,jp,y);
    if (jj == 0 and tmp >= 1.0e-6) { total += tmp; }
    snl += tmp * sint(bn,bex,rdbex,sex,betan,betan.size(),al,wt,tbart,nbx);

  }
  return snl;

}
