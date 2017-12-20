#include "formf.h" 
#include "smallFuncs.h" 


auto fccLatticeFactors( const double& twothd, const double& twopis, 
  int lat, std::vector<double>& b, int ifl, double w, int nw,
  double t2, double c1, double wint, double ulim, double a ){
    
  // compute lattice factors for fcc lattices
  double phi = ulim / twopis, tau, tsq;
  int i1m = a*sqrt(phi), k = 0;
  i1m = 15;
  for ( int i1 = -i1m; i1 <= i1m; ++i1 ){
    for ( int i2 = -i1m; i2 <= i1m; ++i2 ){ 
      for ( int i3 = -i1m; i3 <= i1m; ++i3 ){
        tsq = taufcc( i1, i2, i3, c1, twothd, twopis );
        if (tsq > 0 and tsq <= ulim) {
          k += 1;
          if ((2*k) > nw) std::cout << "storage exceeded" << std::endl; 
          tau = sqrt(tsq);
          b[ifl+2*k-2-1] = tsq;
          b[ifl+2*k-1-1] = ( exp(-tsq*t2*wint) / tau ) * formf(lat,i1,i2,i3);
        }
      }
    }
  }
  return k - 1;
}


