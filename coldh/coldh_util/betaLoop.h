#include <iostream>
#include <vector>
#include <cmath>
#include "betaLoop_util/jprimeLoop.h"
#include "betaLoop_util/bt.h"

auto betaLoop( const std::vector<double>& betan, 
  const std::vector<double>& rdbex, const std::vector<double>& bex,
  const std::vector<double>& sex, const double& alphaVal, const double& wt, 
  const double& tbart, const double& x, const double& y, const double& swe, 
  const double& swo, int itemp, int nbx, int a, int ncold, bool free, 
  std::vector<std::vector<std::vector<double>>>& sym_sab, 
  std::vector<std::vector<std::vector<double>>>& sym_sab_2 ){
  //--loop over all beta values
  //    results for positive beta go into ssp
  //    results for negative beta go into sym_sab
  double pj;
  int jjmax = 2 * betan.size() - 1;
  int k, jp, j;
  double add, betap, be, snlg, snlk, sn;
  double pi = 3.14159265358979;
  for ( auto jj = 0; jj < jjmax; ++jj ){
    if ( jj < betan.size() ){
      k = betan.size() - jj;
      be = betan[k-1];
    } else {
      k = jj - betan.size() + 2;
      be = betan[k-1];
    }
    if ( jj < betan.size()-1 ){
      be = -be;
    }
    double total = 0;
    sn = 0;
    //--loop over all oscillators
    // para-h2: j=0,2,....; ortho-h2: j=1,3,....
    // ortho-d2: j=0,2,....; para-d2: j=1,3,....
    int ipo = ncold == 1 or ncold == 4 ? 2 : 1;

    //int jt1 = 2 * 3; // THIS IS REALLY WEIRD origialy 2 * jterm but jterm
                       // always equals 3? Please check out
    int jt1 = ncold == 1 or ncold == 4 ? 7 : 6;

    for ( auto l = ipo; l < jt1; l = l + 2 ){
      j = l - 1;

      bt(j,pj,x);
      // Get even
      snlg = jPrime( total, j, be, x, swe, pj, jj, bex, rdbex, sex, betan, alphaVal, wt, tbart, y, nbx, false, free );
      // Get odd
      snlk = jPrime( total, j, be, x, swo, pj, jj, bex, rdbex, sex, betan, alphaVal, wt, tbart, y, nbx, true,  free );

      //--continue the j loop
      sn = sn + snlg + snlk;
    }

    //--continue the beta loop
      if (jj < betan.size()) sym_sab[a][k-1][itemp] = sn;
      if (jj >= betan.size()-1) sym_sab_2[a][k-1][itemp] = sn;
  }

}
