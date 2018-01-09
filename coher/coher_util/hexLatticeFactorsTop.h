#include <iostream>
#include <vector>
#include "hexLatticeFactorsMid.h"

auto hexLatticeFactorsTop( double a, double tsq, double c1, double c2, 
  int lat, int nw, double tsqx, std::vector<double>& b, int ifl, 
  int i, double wint, double t2, double ulim, 
  int imax, double c, int i1m ){
  double tau, f;
  // compute lattice factors for hexagonal lattices
  double phi=ulim/(4*M_PI*M_PI), w, w1, w2, w3;
  int l1, l2, l3, i2m, i3m, k = 0;

  for ( auto i1 = 1; i1 <= i1m; ++i1 ){

    l1=i1-1;
    i2m=int((l1+sqrt(3*(a*a*phi-l1*l1)))/2);
    i2m=i2m+1;

    hexLatticeFactorsMid( a, tsq, c1, c2, lat, nw, tsqx, b, ifl,
      i, wint, t2, ulim, imax, c, i1, i2m, l1, k );

  } // 1


  imax = k - 1;
  //go to 220
}


