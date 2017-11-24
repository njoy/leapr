#include <iostream>
#include <vector>
#include "../../discre/discre_util/sint.h"

auto jPrimeOdd(double& total, const double& alp, int j, int ifree, 
    const double& be, double& x, double& swo, double& pj, int jj,
    std::vector<double>& bex, std::vector<double>& rdbex, 
    std::vector<double>& sex, std::vector<double>& betan, double& al,
    double& wt, double& tbart, int maxbb, int b ){
  //--sum over the odd values of j-prime
  double snlk=0;
  double betap;
  double pi = 3.14159265358979;
  double add;
  int jp;
  for ( auto lp = 2; lp < 10; lp = lp + 2 ){
    jp=lp-1;
    betap=(-j*(j+1)+jp*(jp+1))*x/2;
    auto tmp=(2*jp+1)*pj*swo*4;//*sumh(j,jp,y);
    std::cout << tmp << std::endl;
    if (jj == 1 and tmp >= 1.0e-6) {
      total += tmp;
    }
    double bn = be + betap;
    if (ifree == 1) {
      double ex=-(alp-abs(bn))*(alp-abs(bn))/(4*alp);
      if (bn > 0.0) ex=ex-bn;
      add = exp(ex) / sqrt(4*pi*alp);
    } else {
      add = sint(bn,bex,rdbex,sex,betan,b,al,wt,tbart,maxbb);
    }
    snlk = snlk + tmp * add;
  }
  return snlk;

}
