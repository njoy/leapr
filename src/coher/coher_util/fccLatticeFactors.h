#include "coher/coher_util/formf.h" 
#include "generalTools/constants.h"

// FCC
double taufcc( int l1, int l2, int l3, double c1 ){
  return c1 * ( l1*l1 + l2*l2 + l3*l3 + 0.6666667*( l1*l2 + l1*l3 - l2*l3 ) ) * 
	 4 * M_PI * M_PI;
}

auto fccLatticeFactors( int lat, double a, double maxTauSq, double massScatterer,
                        std::vector<double>& b  ){
  double c1=3/(a*a);
    
  // compute lattice factors for fcc lattices
  double phi = maxTauSq / (4*M_PI*M_PI), tau, tsq;
  int i1m = a*sqrt(phi), k = 0;
  i1m = 15;
  for ( int i1 = -i1m; i1 <= i1m; ++i1 ){
    for ( int i2 = -i1m; i2 <= i1m; ++i2 ){ 
      for ( int i3 = -i1m; i3 <= i1m; ++i3 ){
        tsq = taufcc( i1, i2, i3, c1 );
        if (tsq > 0 and tsq <= maxTauSq) {
          k += 1;
          if ((2*k) > int(b.size())) std::cout << "storage exceeded" << std::endl; 
          tau = sqrt(tsq);
          b[2*k-2] = tsq;
          //b[2*k-1] = ( exp(-tsq*t2*wint) / tau ) * formf(lat,i1,i2,i3);
          b[2*k-1] = ( 1.0                 / tau ) * formf(lat,i1,i2,i3);
        }
      }
    }
  }
  return k - 1;
  std::cout << massScatterer << std::endl;
}


