#include <iostream>
#include <vector>
#include "formf.h" 
#include "smallFuncs.h" 


auto bccLatticeFactors( const double& ulim, const double& twopis, 
  int& k, std::vector<double>& b,
  const int& ifl, const double& wint, int& imax, const double& t2, 
  const int& lat, const double& a, const double& c1 ){
  // compute lattice factors for bcc lattices
  //215 continue
  int i1m, i2m, i3m;
  double tau, phi, w, f, tsq;
  phi = ulim/twopis;
  i1m = int(a*sqrt(phi));
  i1m = 15;
  k = 0;
  for ( auto i1 = -i1m; i1 <= i1m; ++i1 ){
    i2m = i1m;
    for (auto i2 = -i2m; i2 <= i2m; ++i2 ){
      i3m = i1m;
      for ( auto i3 = -i3m; i3 <= i3m; ++i3 ){
        tsq = taubcc(i1,i2,i3,c1,twopis);

        if (tsq > 0 and tsq <= ulim) {
          tau = sqrt(tsq);
          w = exp(-tsq*t2*wint)/tau;
          f = w*formf(lat,i1,i2,i3);
          k = k + 1;
          if ((2*k) > b.size() ) std::cout << "storage exceeded" << std::endl;
          b[ifl+2*k-2-1] = tsq;
          b[ifl+2*k-1-1] = f;
        }
      }
    }
  }
  imax=k-1;
}


