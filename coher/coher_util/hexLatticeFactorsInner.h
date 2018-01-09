#include <iostream>
#include <vector>
#include "hexLatticeFactors_util/hexLatticeFactorsHelper.h"
#include "formf.h"
#include "smallFuncs.h"


auto hexLatticeFactorsInner( double& a, double& c1, double& c2, 
  int& lat, int& nw, double& tsqx, std::vector<double>& b, int& ifl, 
  int& i, double& wint, double& t2, double& ulim,
  int& l1, int& l2, int& i3m, int& k ){

  // compute lattice factors for hexagonal lattices
  
  double w, w1, w2, w3, tsq, tau, f;
  int l3; 

  for ( auto i3 = 1; i3 <= i3m; ++i3 ){

    l3 = i3 - 1;
    w1 = 2;
    if (l1 == l2) w1=1;
    w2 = 2;
    if (l1 == 0 or  l2 == 0) w2 = 1;
    if (l1 == 0 and l2 == 0) w2 = w2/2;
    w3 = 2;
    if (l3 == 0) w3 = 1;

    tsq = tausq(l1,l2,l3,c1,c2);

    if (tsq > 0 and tsq <= ulim) {
      tau = sqrt(tsq);
      w = exp(-tsq*t2*wint)*w1*w2*w3/tau;
      f = w*formf(lat,l1,l2,l3);
      hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, wint, nw, f, i );
    }

    tsq = tausq(l1,-l2,l3,c1,c2);

    if (tsq > 0 and tsq <= ulim) {
      tau = sqrt(tsq);
      w = exp(-tsq*t2*wint)*w1*w2*w3/tau;
      f = w*formf(lat,l1,-l2,l3);
      hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, wint, nw, f, i );
    }

  } // 3

}

