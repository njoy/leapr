#include <iostream>
#include <vector>

auto hexLatticeFactorsHelper( int& k, const double& tsq, 
  const double& tsqx, std::vector<double>& b, const int& ifl,
  const double& wint, const int& nw, const double& f ){

  if (k <= 0 or tsq <= tsqx) {
    k += 1;
    if ((2*k) > nw) std::cout << "ERROR" << std::endl; 

      if ( abs(1.80558e-8 - f) < 1e-11 ){ std::cout << tsq << "    " << b[ifl+2*k-3] << "     " << k << "    " << ifl+2*k-3<<  std::endl; }

    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;
  }

  else {

    int i = 0;
    while ( i < k ){
      i += 1;

      //  if ( abs(1.80558e-8 - f) < 1e-11 ){ std::cout << tsq << "    " << b[ifl+2*i-3] << "     " << i << "    " << ifl+2*i-3<<  std::endl; }
      if ( tsq > b[ifl+2*i-3] and tsq < 1.05 * b[ifl+2*i-3] ) {
        b[ifl+2*i-2] += f;

        if ( abs(1.80558e-8 - f) < 1e-11 ){ std::cout << " OH BABY "<< f  << std::endl; }
        return;
      } 
    }

    k += 1;
    if ((2*k) > nw) std::cout << "storage exceeded" << std::endl;
    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;

  } 
} 


