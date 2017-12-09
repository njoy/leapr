#include <iostream>
#include <vector>
#include <cmath>
#include "betaLoop_util/jprimeLoop.h"
#include "betaLoop_util/bt.h"

auto betaLoop( const std::vector<double>& betan, 
  const std::vector<double>& rdbex, const std::vector<double>& bex,
  const std::vector<double>& sex, const double& alphaVal, const double& wt, 
  const double& tbart, const double& x, const double& y, const double& evenSumConst, 
  const double& oddSumConst, int itemp, int nbx, int a, int ncold, bool free, 
  std::vector<std::vector<std::vector<double>>>& sym_sab, 
  std::vector<std::vector<std::vector<double>>>& sym_sab_2 ){

  //--loop over all beta values
  //    results for positive beta go into ssp
  //    results for negative beta go into sym_sab
  int jjmax = 2 * betan.size() - 1, k, jp, j;
  double pj, add, betap, be, snlg, snlk, sn, total;

  for ( auto jj = 0; jj < jjmax; ++jj ){

    k  = jj < betan.size() ?      betan.size() - jj : jj - betan.size() + 2;
    be = jj < betan.size() - 1 ? -betan[k-1]        : betan[k-1];

    //--loop over all oscillators
    // para-h2 : j=0,2,....; ortho-h2: j=1,3,....
    // ortho-d2: j=0,2,....; para-d2 : j=1,3,....
    int start = ncold == 1 or ncold == 4 ? 1 : 0;
    int end   = ncold == 1 or ncold == 4 ? 6 : 5;

    total = 0; sn = 0;

    // This is starting the first sum in Eq. 567-568. For Spara we sum over 
    // even values (Eq. 567), while for Sortho its over odd values (Eq. 568).
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------
    // ALERT ---- One problem with this is that para-h2 and ortho-d2 are being
    // treated the same, just like they are when being multiplied by sk in
    // coldh.h
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------

    for ( auto j = start; j < end; j = j + 2 ){

      bt(j,pj,x);

      // Perform inner sum over even values (inner sum in Eq. 567-568 that is
      // multiplied by A)
      snlg = jPrime( total, j, be, x, evenSumConst, pj, jj, bex, rdbex, sex, betan, 
        alphaVal, wt, tbart, y, nbx, false, free );

      // Perform inner sum over odd values (inner sum in Eq. 567-568 that is
      // multiplied by B)
      snlk = jPrime( total, j, be, x, oddSumConst, pj, jj, bex, rdbex, sex, betan, 
        alphaVal, wt, tbart, y, nbx, true,  free );

      //--continue the j loop
      sn = sn + snlg + snlk;
    }

    //--continue the beta loop
      if (jj < betan.size())    sym_sab[a][k-1][itemp]   = sn;
      if (jj >= betan.size()-1) sym_sab_2[a][k-1][itemp] = sn;
  }

}
