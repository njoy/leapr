#include "simple/continuous/convolution.h"
#include <iostream>


double fact(int n){
  if (n == 0){ return 1; }
  return n*fact(n-1);
}


void print(std::vector<double> x){
  std::cout << std::endl;
  for (auto& y : x){ std::cout << y << " "; }
  std::cout << std::endl;
}
void print(double x){
  std::cout << std::endl;
  std::cout << x << std::endl;
}





template <typename V, typename F>
auto sumOverTn( const V& alphas, const V& betas, const V& T1, 
  const F& lambda_s ){
  alphas = {2,4};
  betas = {0,1,2,3,4};
  T1 = {1,2,3,2,1};
  lambda_s = 0.5;
  int nAlpha = alphas.size();
  int nBeta  = betas.size();


  nAlpha = alphas.size();
  nBeta  = betas.size();
  V Tn(nBeta);
  std::copy( T1.begin(), T1.begin() + nBeta, Tn.begin() );



  F factorialTerm = 1.0;
  V sab(nAlpha*nBeta,0.0);
  // Perform sum in Eq. 523
  for ( int n = 0; n < 2; ++n ){
    for ( int a = 0; a < nAlpha; ++a ){
      F alpha = alphas[a];
      F alphaTerm = exp(-alpha*lambda_s)*(1.0/fact(n))*std::pow(alpha*lambda_s,n);
      if (n == 0 and betas[0] == 0){
        sab[a*nBeta+0] = alphaTerm;
        continue;
      }
      for ( int b = 0; b < nBeta; ++b ){
        sab[a*nBeta+b] += Tn[b]*alphaTerm;
      } // beta loop
    } // alpha loop
    if (n == 0){ continue; }
    Tn = convolve(betas,T1,Tn);
  } // summation loop

  //print(sab[0*nBeta+0]);
  //print(sab[0*nBeta+1]);
  print(exp(-1)*3);
  print(sab[0*nBeta+2]);





}


