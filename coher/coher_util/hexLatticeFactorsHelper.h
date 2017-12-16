#include <iostream>
#include <vector>
#include "smallFuncs.h"
#include "formf.h"


auto hexLatticeFactorsHelper( int& k, const double& tsq, 
  const double& tsqx, std::vector<double>& b, const int& ifl,
  const double& wint, const int& nw, const double& eps, 
  double f ){

  double tau;
  if (k <= 0 or tsq <= tsqx) {
    k=k+1;
    if ((2*k) > nw) std::cout << "ERROR" << std::endl; 
    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;
  }
  else {
    int i = 0;
    bool done = false;
    while (i < k and not done ){
      i += 1;
      if (tsq >= b[ifl+2*i-3] and tsq < (1+eps)*b[ifl+2*i-3]) {
        b[ifl+2*i-2] += f;
        done = true;
      } 
    }
    if (not done) {
      k += 1;
      if ((2*k) > nw) std::cout << "storage exceeded" << std::endl;
      b[ifl+2*k-3] = tsq;
      b[ifl+2*k-2] = f;
    }
  } 
} 


