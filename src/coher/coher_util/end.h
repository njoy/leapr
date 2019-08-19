#include "generalTools/swap.h"


auto sortLatticeFactors(std::vector<double>& b, int& k, double maxTauSq, int imax){
  for ( auto i = 0; i < imax; ++i ){
    for ( auto j = i+1; j < k; ++j ){
      if (b[2*j] < b[2*i]) {
	swap(b[2*i  ],b[2*j  ]);
	swap(b[2*i+1],b[2*j+1]);
      }
    }
  }
  k++;
  b[2*k-2] = maxTauSq;
  b[2*k-1] = b[2*k-3];
  //nw = 2 * k;
}


auto end( std::vector<double>& b, int& k, double econ, 
  double toler, double scon, double maxTauSq, int imax ){

  sortLatticeFactors( b, k, maxTauSq, imax );

  // convert to practical units and combine duplicate bragg edges.
  double bel = -1, be, bs, recon = 1.0/econ;
  int nbe, j = 0;

  for ( auto i = 1; i <= k; ++i ){
    be = b[2*i-2] * recon;
    bs = b[2*i-1] * scon;
    if (be-bel < toler) {
      b[2*j-1] = b[2*j-1] + bs;
    }
    else {
      j = j+1;
      b[2*j-2] = be;
      b[2*j-1] = bs;
      bel=be;
    }
  }

  nbe = j;
  // only invoke this for iel = 2, don't forget you're debugging weird 
  // stuff in formf right now
  return nbe;
}


