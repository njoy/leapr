#include <iostream>
#include <vector>
#include "../../discre/discre_util/sint.h"
#include "sumh.h"

auto jPrimeOdd(double& total, int j, int ifree, 
  const double& be, double& x, double& swo, double& pj, int jj,
  std::vector<double>& bex, std::vector<double>& rdbex, 
  std::vector<double>& sex, std::vector<double>& betan, double& al,
  double& wt, double& tbart, int b, double& y, int nbx ){

  //--sum over the odd values of j-prime
  
  double betap, add, snlk = 0, pi = 3.14159265358979, tmp, bn;
  int jp;

  for ( auto jp = 1; jp < 10; jp = jp + 2 ){
    betap = ( -j * (j+1) + jp * (jp+1) ) * x * 0.5;
    tmp = (2*jp+1) * pj * swo * 4 * sumh(j,jp,y);
    if (jj == 0 and tmp >= 1.0e-6) {
      total += tmp;
    }
    bn = be + betap;
    snlk = snlk + tmp * sint(bn,bex,rdbex,sex,betan,b,al,wt,tbart,nbx);

  }
  return snlk;

}
auto jPrimeEven(double& total, int j, int ifree, 
  const double& be, double& x, double& swe, double& pj, int jj,
  std::vector<double>& bex, std::vector<double>& rdbex, 
  std::vector<double>& sex, std::vector<double>& betan, double& al,
  double& wt, double& tbart, int b, double& y, int nbx ){


  //--sum over even values of j-prime
  double betap, add, snlk = 0, pi = 3.14159265358979, tmp, bn;
  int jp;

  double snlg = 0.0;
  for ( auto jp = 0; jp < 9; jp = jp + 2 ){
    betap = ( -j * (j+1) + jp * (jp+1) ) * x * 0.5;
    tmp = (2*jp+1) * pj * swe * 4 * sumh(j,jp,y);
    if (jj == 0 and tmp >= 1.0e-6) {
      total += tmp;
    }
    bn = be + betap;
    snlg = snlg + tmp * sint(bn,bex,rdbex,sex,betan,b,al,wt,tbart,nbx);
  }
  return snlg;
} 
