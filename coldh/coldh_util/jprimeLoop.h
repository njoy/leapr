#include <iostream>
#include <vector>
#include "../../discre/discre_util/sint.h"
#include "sumh.h"

auto jPrimeOdd(double& total, const double& alp, int j, int ifree, 
    const double& be, double& x, double& swo, double& pj, int jj,
    std::vector<double>& bex, std::vector<double>& rdbex, 
    std::vector<double>& sex, std::vector<double>& betan, double& al,
    double& wt, double& tbart, int maxbb, int b, double& y, int nbx ){
  //--sum over the odd values of j-prime
  double betap, add, snlk = 0, pi = 3.14159265358979;
  int jp;
  for ( auto jp = 1; jp <= 9; jp = jp + 2 ){
    betap=(-j*(j+1)+jp*(jp+1))*x/2;
    auto tmp=(2*jp+1)*pj*swo*4*sumh(j,jp,y);
    if (jj == 0 and tmp >= 1.0e-6) {
      total += tmp;
    }
    double bn = be + betap;
    // sooo this ifree thing is never not equal to zero
    if (ifree == 1) {
      double ex=-(alp-abs(bn))*(alp-abs(bn))/(4*alp);
      if (bn > 0.0) { ex = ex - bn; }
      add = exp(ex) / sqrt(4*pi*alp);
    } else {
      add = sint(bn,bex,rdbex,sex,betan,b,al,wt,tbart,nbx);
    }
    snlk = snlk + tmp * add;
  }
  return snlk;

}
