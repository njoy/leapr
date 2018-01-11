#include <iostream>
#include <vector>

auto hexLatticeFactorsHelper( int& k, const double& tsq, 
  const double& tsqx, std::vector<double>& b, const int& ifl,
  const double& wint, const int& nw, const double& f, int& i ){
 
  if (k <= 0 or tsq <= tsqx) {
    k += 1;

    if ((2*k) > nw) std::cout << "ERROR" << std::endl; 

    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;
  }

  else {
    i = 0;

    while ( i < k ){
      i += 1;

      if ( tsq >= b[ifl+2*i-3] and tsq <= 1.05 * b[ifl+2*i-3] ) {
        b[ifl+2*i-2] += f;
	return;
      } // if
    } // while

    k += 1;
    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;

  } // else
}
