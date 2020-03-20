#include "coher/coher_util/formf.h" 
#include "generalTools/constants.h"

// BCC
template <typename Float>
Float taubcc( int l1, int l2, int l3, Float c1 ){
  return c1 * ( l1*l1 + l2*l2 + l3*l3 + l1*l2 + l2*l3 + l1*l3 ) *
          4 * M_PI * M_PI;
}

template <typename Float>
auto bccLatticeFactors( const Float& maxTauSq, std::vector<Float>& b, 
  const int& lat, const Float& a, const Float& massScatterer ){

  double c1 = 2.0/(a*a), tsq;
  // compute lattice factors for bcc lattices
  int im;
  unsigned int k = 0;
  im = 15; // = int(a*sqrt(maxTauSq/(4*M_PI*M_PI)));

  for ( auto i1 = -im; i1 <= im; ++i1 ){
    for ( auto i2 = -im; i2 <= im; ++i2 ){
      for ( auto i3 = -im; i3 <= im; ++i3 ){
        tsq = taubcc(i1,i2,i3,c1);
        if (tsq > 0 and tsq <= maxTauSq) {
          k += 1;
          if ((2*k) > b.size() ) std::cout << "storage exceeded (BCC)" << std::endl;
          b[2*k-2] = tsq;
          //b[2*k-1] = exp(-tsq*t2*wint)*formf(lat,i1,i2,i3)/sqrt(tsq);
          b[2*k-1] =                     formf(lat,i1,i2,i3)/sqrt(tsq);
        }
      } // i3
    } // i2 
  } // i1

  return k-1; // this is imax
  std::cout << massScatterer << std::endl;
}

