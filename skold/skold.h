#include <iostream>
#include <vector>
#include "../coldh/coldh_util/terpk.h"
#include "skold_util/terp1.h"

auto skold( double cfrac, int itemp, double tev,
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<double>& skappa, int ntempr, double awr, int lat, int nka, 
  double dka, double scaling,
  std::vector<std::vector<std::vector<double>>>& symSab ){
  /* use skold approximation to add in the effects
   * of intermolecular coherence.
   */
  int i, j, k, kk, nal, ibeta;
  double sc, sk, ap, be, ss, s1, s2, waven;
  std::vector<double> scoh ( 1000, 0.0 );

  double amassn = 1.008664904,    amu = 1.6605402e-24, bk = 8.617385e-5,
         hbar   = 1.05457266e-27, ev  = 1.60217733e-12;

  // apply the skold approximation
  sc = 1;
  if (lat == 1) sc = 0.0253 / tev;
  for ( auto b = 0; b < beta.size(); ++b ){
    for ( auto a = 0; a < alpha.size(); ++a ){
      waven = 1.0e-8 * sqrt(2*awr*amassn*amu*tev*ev*alpha[a]*scaling)/hbar;
      sk = terpk(skappa,dka,waven);
      ap = alpha[a] / sk;
      for ( auto a2 = 0; a2 < alpha.size(); ++a2 ){
        kk = a2;
        if (ap < alpha[a2]){ break; }
      }
      if (kk == 0) kk = 1;

      terp1( alpha[kk-1], symSab[kk-1][b][itemp], alpha[kk], 
             symSab[kk][b][itemp], ap,s coh[a], 5 );
      scoh[a] *= sk;
    }
    for ( auto a = 0; a < alpha.size(); ++a ){
      symSab[a][b][itemp] = (1-cfrac)*symSab[a][b][itemp]+cfrac*scoh[a];
    }
  }
}

