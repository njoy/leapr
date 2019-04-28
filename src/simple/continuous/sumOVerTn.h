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



template <typename V>
void addPadding( V& betas, V& T1, V& T2, int numZeros ){

  // Find a good representative value of delta beta so that I can expand it. 
  // For now, just choose the end delta
  V pad(numZeros,0.0);
  auto Delta = betas[1]-betas[0];
  for (int i = 0; i < numZeros; ++i){
    pad[i] = betas[betas.size()-1]+(i+1)*Delta;
  }
  betas.insert( betas.end(), pad.begin(), pad.end() );
  std::reverse(pad.begin(),pad.end());
  for (int i = 0; i < numZeros; ++i){ pad[i] *= -1; }

  betas.insert( betas.begin(), pad.begin(), pad.end() );

  V zeros (numZeros,0.0);
  T1.insert( T1.end(), zeros.begin(), zeros.end() );
  T1.insert( T1.begin(), zeros.begin(), zeros.end() );
  T2.insert( T2.end(), zeros.begin(), zeros.end() );
  T2.insert( T2.begin(), zeros.begin(), zeros.end() );
}




template <typename V, typename F>
V sumOverTn( const V& alphas, V betas, V T1, const F& lambda_s,
  int N = 0){

  int nAlpha = alphas.size();
  int nBeta  = betas.size();

  V Tn(nBeta);
  std::copy( T1.begin(), T1.begin() + nBeta, Tn.begin() );

  V beta_full = reflect(betas), T1_full = reflect(T1);
  for ( int i = 0; i < int(beta_full.size()/2); ++i ){ 
    beta_full[i] *= -1; 
    T1_full[i] *= exp(-beta_full[i]); 
  }

  V beta_original = beta_full;
  V T2_full(beta_full.size());
  std::copy(T1_full.begin(), T1_full.begin() + beta_full.size(), T2_full.begin());

  // This is going to add terms to my beta ends, and add zeros to my T vecs
  int numZeros = 1;
  addPadding(beta_full,T1_full,T2_full,numZeros);

  F factorialTerm = 1.0;
  V sab(nAlpha*nBeta,0.0);
  // Perform sum in Eq. 523
  for ( int n = 0; n < N; ++n ){
    for ( int a = 0; a < nAlpha; ++a ){
      F alpha = alphas[a];
      F alphaTerm = exp(-alpha*lambda_s)*(1.0/fact(n))*std::pow(alpha*lambda_s,n);

      // Accounting for our T0(b) = delta(b) term
      if (n == 0){
        if (betas[0] == 0){ sab[a*nBeta+0] = alphaTerm; } 
        continue;
      }

      // now i need to add in all the relevant terms into my S(a,b)
      int halfway = T2_full.size()*0.5;
      for ( int b = 0; b < nBeta; ++b ){
        sab[a*nBeta+b] += Tn[b]*alphaTerm;
        //sab[a*nBeta+b] += T2_full[halfway+b]*alphaTerm;
      } // beta loop

    } // alpha loop
    if (n == 0){ continue; }
    Tn = convolve(betas,T1,Tn);
    T2_full = convolutionWithPadding(beta_original,beta_full,T1_full,T2_full);

  } // summation loop
  return sab;

}



















