#include "simple/continuous/convolution.h"
#include <iostream>


double fact(int n){ return n == 0 ? 1 : n*fact(n-1); }

void print(std::vector<double> x){
  std::cout << std::endl;
  for (auto& y : x){ std::cout << y << " "; }
  std::cout << std::endl;
}

template<typename T>
void print(T x){
  std::cout << std::endl << x << std::endl;
}
template<typename T>
void print(T x,T y){
  std::cout << std::endl << x  << "    " << y<< std::endl;
}




template <typename V, typename F>
V sumOverTn( const V& alphas, V betas, V T1, const F& lambda_s,
  int N = 0){

  int nAlpha = alphas.size();
  int nBeta  = betas.size();

  V Tn(nBeta);
  std::copy( T1.begin(), T1.begin() + nBeta, Tn.begin() );



  F factorialTerm = 1.0;
  V sab(nAlpha*nBeta,0.0);
  // Perform sum in Eq. 523
  for ( int n = 0; n < N; ++n ){
    for ( int a = 0; a < nAlpha; ++a ){
      F alpha = alphas[a];
      F alphaTerm = exp(-alpha*lambda_s)*(1.0/fact(n))*std::pow(alpha*lambda_s,n);
      //if ( a == 0){ print(1.0*n,alphaTerm); }
      // Accounting for our T0(b) = delta(b) term
      if (n == 0){
        if (betas[0] == 0){ sab[a*nBeta+0] = alphaTerm; } 
        continue;
      }
      for ( int b = 0; b < nBeta; ++b ){
        sab[a*nBeta+b] += Tn[b]*alphaTerm;
      } // beta loop
    } // alpha loop
    if (n == 0){ continue; }
    Tn = convolve(betas,T1,Tn);
  } // summation loop

  return sab;

}


