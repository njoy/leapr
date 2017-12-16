#include <iostream>
#include <vector>
#include "smallFuncs.h"
#include "formf.h"


auto hexagonalLatticeFactorsHelper( const int& lat, const int& l1, const int& l2, 
  const int& l3, double& w, int& k, const double& tsq, const double& tsqx, 
  std::vector<double>& b, const int& ifl, int& i,
  const double& ulim, const int& t2, const double& w1, const double& w2, 
  const double& w3, const double& wint, const int& nw, const double& eps ){

  double f, tau;
  tau=sqrt(tsq);
  w=exp(-tsq*t2*wint)*w1*w2*w3/tau;
  f=w*formf(lat,l1,l2,l3);
  if (k <= 0 or tsq <= tsqx) {
    k=k+1;
    if ((2*k) > nw) std::cout << "ERROR" << std::endl; 
    b[ifl+2*k-2-1]=tsq;
    b[ifl+2*k-1-1]=f;
  }
  else {
    i=0;
    int idone=0;
    while (i < k and idone == 0){
      i=i+1;
      if (tsq >= b[ifl+2*i-2-1] and tsq < (1+eps)*b[ifl+2*i-2-1]) {
        b[ifl+2*i-1-1]=b[ifl+2*i-1-1]+f;
        idone=1;
      } 
    }
    if (idone == 0) {
      k=k+1;
      if ((2*k) > nw) std::cout << "storage exceeded" << std::endl;
      b[ifl+2*k-2-1]=tsq;
      b[ifl+2*k-1-1]=f;
    }
  } 
} 


