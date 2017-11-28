#include <iostream>
#include <vector>
#include "../../discre/discre_util/sint.h"
#include "sumh.h"

auto jPrime( double& total, int j, int ifree, 
  const double& be, const double& x, const double& sw, const double& pj, 
  int jj, const std::vector<double>& bex, const std::vector<double>& rdbex, 
  const std::vector<double>& sex, const std::vector<double>& betan, 
  const double& al, const double& wt, const double& tbart, int b, 
  const double& y, int nbx, bool odd ){

  //--sum over the odd or even values of j-prime
  
  double betap, add, snl = 0, pi = 3.14159265358979, tmp, bn;
  int jp, start, end;
  if ( odd == true ){ start = 1; end = 10; }
  else              { start = 0; end = 9;  }

  for ( auto jp = start; jp < end; jp = jp + 2 ){
    betap = ( -j * (j+1) + jp * (jp+1) ) * x * 0.5;
    tmp = (2*jp+1) * pj * sw * 4 * sumh(j,jp,y);
    if (jj == 0 and tmp >= 1.0e-6) { total += tmp; }
    bn = be + betap;
    snl += tmp * sint(bn,bex,rdbex,sex,betan,b,al,wt,tbart,nbx);

  }
  return snl;

}
