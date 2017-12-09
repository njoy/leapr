#include <iostream>
#include <vector>
#include "jprimeLoop_util/sumh.h"
#include "../../../discre/discre_util/sint.h"

auto jPrime( double& total, int j, const double& be, const double& x, 
  const double& sw, const double& pj, int jj, const std::vector<double>& bex, 
  const std::vector<double>& rdbex, const std::vector<double>& sex, 
  const std::vector<double>& betan, const double& al, const double& wt, 
  const double& tbart, const double& y, int nbx, bool odd, bool free  ){

  //--sum over the odd or even values of j-prime
  
  double betap, add, snl = 0, pi = 3.14159265358979, tmp, bn, ex;
  int jp, start, end;

  // Select start and end values for loop. If we want odd values, 1-->9,
  // if even 2-->8.
  if ( odd == true ){ start = 1; end = 10; }
  else              { start = 0; end = 9;  }

  for ( auto jp = start; jp < end; jp = jp + 2 ){

    bn = be + ( -j*(j+1) + jp*(jp+1) ) * x * 0.5;
    tmp = (2*jp+1) * pj * sw * 4 * sumh(j,jp,y);

    if (jj == 0 and tmp >= 1.0e-6) { total += tmp; }

    if ( free ) {
      ex = -std::pow(al*wt-abs(bn),2)/(4*al*wt);
      if ( bn > 0.0 ){ ex -= bn; }
      add = exp(ex)/sqrt(4*pi*al*wt);
    }
    else{
      add = sint(bn,bex,rdbex,sex,betan,betan.size()-1,al,wt,tbart,nbx);
    }

    snl += tmp * add;

  }
  return snl;

}
