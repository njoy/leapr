#include "generalTools/tools.h"


template <typename Range, typename Float>
auto sortLatticeFactors(Range& b, int& k, const Float& maxTauSq, int imax){
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
}


template <typename Range, typename Float>
auto end( Range& b, int& k, const Float& econ, const Float& toler, 
  const Float& scon, const Float& maxTauSq, int imax ){

  sortLatticeFactors( b, k, maxTauSq, imax );

  // convert to practical units and combine duplicate bragg edges.
  Float bel = -1, be, bs, recon = 1.0/econ;
  int j = 0;

  for ( int i = 1; i <= k; ++i ){
    be = b[2*i-2] * recon;
    bs = b[2*i-1] * scon;
    if (be-bel < toler) {
      b[2*j-1] += bs;
      continue;
    }
    ++j;
    b[2*j-2] = be;
    b[2*j-1] = bs;
    bel = be;
  }
  return j; // j = nbe
}


