#include <iostream>
#include <vector>
#include "formf.h" 
#include "smallFuncs.h" 


auto bccLatticeFactors( double& phi, double& ulim, double& twopis, int& i1m,
  int& k, int& i2m, int& i3m, double& tau, double& w, std::vector<double>& b,
  int ifl, double& tsq, double& f, double& wint, int imax, double& t2, 
  int lat, int nw, double& a, double& c1 ){
  // compute lattice factors for bcc lattices
  //215 continue
  phi = ulim/twopis;
  i1m = int(a*sqrt(phi));
  i1m = 15;
  k = 0;
  for ( auto i1 = -i1m; i1 < i1m; ++i1 ){
    i2m = i1m;
    for (auto i2 = -i2m; i2 < i2m; ++i2 ){
      i3m = i1m;
      for ( auto i3 = -i3m; i3 < i3m; ++i3 ){
        tsq = taubcc(i1,i2,i3,c1,twopis);
        if (tsq > 0 and tsq <= ulim) {
          tau = sqrt(tsq);
          w = exp(-tsq*t2*wint)/tau;
          f = w*formf(lat,i1,i2,i3);
          k = k + 1;
          if ((2*k) > nw) std::cout << "storage exceeded" << std::endl;
          b[ifl+2*k-2-1] = tsq;
          b[ifl+2*k-1-1] = f;
        }
      }
    }
  }
  imax=k-1;
}


