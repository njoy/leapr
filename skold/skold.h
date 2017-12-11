#include <iostream>
#include <vector>
#include "../coldh/coldh_util/terpk.h"
#include "skold_util/terp1.h"

auto skold( double cfrac, int itemp, double temp, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<double>& skappa,
  int ntempr, double awr, int lat, int nka, double dka, double scaling,
  std::vector<std::vector<std::vector<double>>>& symSab ){
  /* use skold approximation to add in the effects
   * of intermolecular coherence.
   */
  int i, j, k, kk, nal, ibeta;
  double tev, sc, amass, al, sk, ap, be, ss, s1, s2;
  double sum0, sum1, ff1l, ff2l, bel, ff1, ff2, waven;
  std::vector<double> scoh ( 1000, 0.0 );

  double amassn = 1.008664904; 
  double amu = 1.6605402e-24;
  double bk = 8.617385e-5;
  double hbar = 1.05457266e-27;
  double ev = 1.60217733e-12;

  
  // apply the skold approximation
  tev = bk * abs(temp);
  sc = 1;
  if (lat == 1) sc = 0.0253 / tev;
  amass = awr * amassn * amu;
  for ( auto b = 0; b < beta.size(); ++b ){
    for ( auto a = 0; a < alpha.size(); ++a ){
      al = alpha[a] * scaling;
      waven = 1.0e-8 * sqrt(2*amass*tev*ev*al)/hbar;
      sk = terpk(skappa,dka,waven);
      ap = alpha[a] / sk;
      for ( auto a2 = 0; a2 < alpha.size(); ++a2 ){
        kk = a2;
        if (ap < alpha[a2]) break;
      }
      if (kk == 0) kk = 1;
      terp1(alpha[kk-1],symSab[kk-1][b][itemp],alpha[kk],symSab[kk][b][itemp],ap,scoh[a],5);
      scoh[a] = scoh[a] * sk;
    }
    for ( auto a = 0; a < alpha.size(); ++a ){
      symSab[a][b][itemp] = (1-cfrac)*symSab[a][b][itemp]+cfrac*scoh[a];
    }
  }
}

