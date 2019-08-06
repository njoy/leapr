#include <iostream>
#include <vector>
#include <cmath>
#include "betaLoop_util/jprimeLoop.h"
#include "betaLoop_util/bt.h"

template <typename Float, typename Range>
auto betaLoop( const Range& betan, const Range& rdbex, const Range& bex, 
  const Range& sex, const Float& alphaWgt, const Float& tbart, 
  const Float& x, const Float& y, const Float& evenSumConst, 
  const Float& oddSumConst, int nbx, int a, int ncold, bool free, 
  Range& sym_sab_1, Range& sym_sab_2 ){
  
  //--loop over all beta values
  //    results for positive beta go into ssp
  //    results for negative beta go into sym_sab
  unsigned int jjmax = 2 * betan.size() - 1; 
  int k;
  Float pj, be, snlg, snlk, sn;

  for ( size_t jj = 0; jj < jjmax; ++jj ){

    k  = jj < betan.size()     ?  betan.size() - jj : jj - betan.size() + 2;
    be = jj < betan.size() - 1 ? -betan[k-1]        : betan[k-1];

    //--loop over all oscillators
    // para-h2 : j=0,2,....; ortho-h2: j=1,3,....
    // ortho-d2: j=0,2,....; para-d2 : j=1,3,....
    int start = ncold == 1 or ncold == 4 ? 1 : 0;
    int end   = ncold == 1 or ncold == 4 ? 6 : 5;

    sn = 0;

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

      /* This calcualtes Pj, the statistic weighting factor for Eq. 567 - 568.
       * 
       *                  ( 2j + 1 ) * exp( -j*(j+1)*x/2 ) 
       * Pj = -----------------------------------------------------------------
       *                    2 * sum (2*k+1)e^( -k*(k+1)*x/2 )
       *
       * where the sum is over even 1-->19 if j is odd (ortho h or para d), or
       * over 0-->18 if j is even (para h or ortho d)
       *
       */

      pj = bt(j,x);

      // Perform inner sum over even values (inner sum in Eq. 567-568 that is
      // multiplied by A)
      snlg = jPrime( j, be, x, evenSumConst, pj, bex, rdbex, sex, betan, 
        alphaWgt, tbart, y, nbx, false, free );

      // Perform inner sum over odd values (inner sum in Eq. 567-568 that is
      // multiplied by B)
      snlk = jPrime( j, be, x, oddSumConst, pj, bex, rdbex, sex, betan, 
        alphaWgt, tbart, y, nbx, true,  free );

      //--continue the j loop
      sn = sn + snlg + snlk;
    }

    //--continue the beta loop
      if (jj <  betan.size())   sym_sab_1[(k-1)+a*betan.size()] = sn;
      if (jj >= betan.size()-1) sym_sab_2[(k-1)+a*betan.size()] = sn;
  } // end for 

}
