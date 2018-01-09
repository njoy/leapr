#include <iostream>
#include <vector>
#include "hexLatticeFactorsInner.h"


auto hexLatticeFactorsMid( double& a, double& tsq, double& c1, double& c2, 
  int& lat, int& nw, double& tsqx, std::vector<double>& b, int& ifl, 
  int& i, double& wint, double& t2, double& ulim, 
  int& imax, double& c, int& i1, int& i2m, int& l1, int& k ){
  double tau, f;
  // compute lattice factors for hexagonal lattices
  double phi=ulim/(4*M_PI*M_PI), w, w1, w2, w3;
  int l2, l3, i3m;

  for ( auto i2 = i1; i2 <= i2m; ++i2 ){

    l2=i2-1;

    double x = phi-c1*(l1*l1+l2*l2-l1*l2);
    i3m = 0;
    if (x > 0 ) i3m = int(c*sqrt(x));
    i3m = i3m + 1;

    hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, i, wint, t2, 
      ulim, l1, l2, i3m, k );


  } // 2
}


