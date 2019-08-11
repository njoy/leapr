#include "generalTools/constants.h"

template <typename Float, typename Range>
auto hexLatticeFactorsHelper( int& k, const Float& tsq, const Float& tsqx, 
  Range& b, const Float& f ){
 
  if (k <= 0 or tsq <= tsqx) {
    k += 1;
    if (2*k > int(b.size())) { throw std::invalid_argument("2k must be <= b size"); } 
    b[2*k-2] = tsq;
    b[2*k-1] = f;
  }
  else {
    for ( int i = 1; i < k+1; ++i ){
      if ( tsq >= b[2*i-2] and tsq <= 1.05 * b[2*i-2] ) {
        b[2*i-1] += f;
	return;
      } // if
    } // while
    k++;
    b[2*k-2] = tsq;
    b[2*k-1] = f;
  } // else
}

