#include <iostream>
#include <vector>
#include "hexLatticeFactors_util/hexLatticeFactorsHelper.h"
#include "coher/coher_util/formf.h" 

// HEXAGONAL
double tausq( int i1, int i2, int i3, double c1, double c2 ){
  /* This is meant to sort of evaluate Eq. 558. The output is  
   * multiplied by an extra factor of 2pi.
   * This is the reciprocal lattice vector length for hexagonal lattice 
   * (e.g. graphite). 
   */
  return (c1 * (i1*i1 + i2*i2 + i1*i2) + (i3*i3*c2)) * 4 * M_PI * M_PI;
}
// According to the HEXSCAT documentation, this actually calcualtes tau*tau
// because (tau/2pi)^2 = (4/3 a*a)(l1*l1 + l2*l2 + l1*l2) + l3*l3/(c*c)


int hexLatticeFactors( double a, double c1, double c2, 
  int lat, double tsqx, std::vector<double>& b,  
  int i, double wint, double t2, double maxTauSq, double c, int i1m ){
  using std::sqrt;

  double tau, f, tsq, tsq2;
  // compute lattice factors for hexagonal lattices
  double phi=maxTauSq/(4*M_PI*M_PI), w, w1, w2, w3;
  //int i1m_b = a*sqrt(phi);
  int i2m, i3m, k = 0;

  for ( auto i1 = 0; i1 < i1m; ++i1 ){
    i2m = int((i1+sqrt(3*(a*a*phi-i1*i1)))/2) + 1;

    for ( auto i2 = i1; i2 < i2m; ++i2 ){
      double x = phi-c1*(i1*i1+i2*i2-i1*i2);

      i3m = (x > 0) ? int(c*sqrt(x)) + 1 : 1;

      for ( auto i3 = 0; i3 < i3m; ++i3 ){
        w1 = (i1 == i2) ? 1 : 2;
        w3 = (i3 == 0) ? 1 : 2;
        w2 = 2;
        if (i1 == 0 or  i2 == 0) w2 = 1;
        if (i1 == 0 and i2 == 0) w2 = w2/2;

        tsq = tausq(i1,i2,i3,c1,c2);

        if (tsq > 0 and tsq <= maxTauSq) {
          tau = sqrt(tsq);   // w1 w2 w3 are weighting factors
          f = exp(-tsq*t2*wint)*w1*w2*w3/tau * formf(lat,i1,i2,i3);
          hexLatticeFactorsHelper( k, tsq, tsqx, b, f, i );
        }

        tsq = tausq(i1,-i2,i3,c1,c2);

        if (tsq > 0 and tsq <= maxTauSq) {
          tau = sqrt(tsq);   // w1 w2 w3 are weighting factors
          f = exp(-tsq*t2*wint)*w1*w2*w3/tau * formf(lat,i1,-i2,i3);
          hexLatticeFactorsHelper( k, tsq, tsqx, b, f, i );
        }

      } // 3

    } // 2

  } // 1

  return k-1; // This is imax

}


