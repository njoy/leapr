#ifndef LEAPR_SIMPLE_CONTINUOUS_CONVOLUTION
#define LEAPR_SIMPLE_CONTINUOUS_CONVOLUTION


#include "simple/generalTools/interpolate.h"
#include "simple/generalTools/integrate.h"
#include <iostream>




template <typename V>
V convolutionWithPadding(V& beta_original, V& beta, V& T1, V& T2){
  /* I'm going to assume that I've been given beta, T1, and T2 such that there
   * is a reasonable amount of zero's padding each side, so that when I convolve
   * T1 and T2 to get T3, I won't have to worry abount dropping terms off the
   * end. This is something you should be worrying about in the upstairs func.
   * I also care about getting b_original (which is the beta grid that describes
   * where the nonzero values of T1 are) so that I can easily tell which values
   * of beta' to traverse across. Note that no matter what values of T2 are 
   * nonzero, the only values of beta' that will result in nonzero contributions
   * are those that have T1(beta') != 0
   */

  // How many zeros did I add on to go from beta_original --> beta ? 
  int numZeros = 0.5*(beta.size() - beta_original.size());
  
  V T3 (beta.size(),0.0);
  for (size_t b = 0; b < beta.size(); ++b){
    V integrand (beta.size(),0.0);
    for (size_t bp = 0; bp < beta_original.size(); ++bp){
      auto T1_piece = interpolate(beta,T1,beta_original[bp]);
      auto T2_piece = interpolate(beta,T2,beta[b]-beta_original[bp]);
      integrand[numZeros+bp] = T1_piece*T2_piece;
    }
    T3[b] = integrate(beta,integrand);
  }
  return T3;

}










template <typename V>
V reflect(V v){
  V fullVec( 2*v.size()-1 );
  std::reverse_copy (v.begin(), v.begin()+v.size(), fullVec.begin());
  std::copy (v.begin(), v.begin()+v.size(), fullVec.begin()+int(fullVec.size()/2));
  return fullVec;
}


template <typename V>
V convolve(V& b, V& y1_0, V y2_0){ 
  V bp = reflect(b);
  V y1 = reflect(y1_0);
  V y2 = reflect(y2_0);
  for ( int i = 0; i < int(bp.size()/2); ++i ){ bp[i] *= -1; }
  for ( int i = 0; i < int(y1.size()/2); ++i ){ y1[i] *= exp(-bp[i]); }
  for ( int i = 0; i < int(y2.size()/2); ++i ){ y2[i] *= exp(-bp[i]); }
  //for (auto& x : y1){ std::cout << x << "  "; }
  //std::cout << std::endl;

  V y3 (b.size()+1,0.0);

  for (size_t i = 0; i < y3.size(); ++i ){
    for (size_t j = 0; j < bp.size()-1; ++j ){
      auto delta = bp[j+1]-bp[j];
      double b_i = (i < b.size()) ?  b[i] : b[i-1]+(b[i-1]-b[i-2]);
      y3[i] += delta/2*(interpolate(bp,y1,bp[j])*interpolate(bp,y2,b_i-bp[j])) + 
               delta/2*(interpolate(bp,y1,bp[j+1])*interpolate(bp,y2,b_i-bp[j+1]));
    }
  }
  b.resize(b.size()+1);
  b[b.size()-1] = b[b.size()-2] + (b[b.size()-2] - b[b.size()-3]);
  y1_0.resize(y1_0.size()+1);
  y1_0[y1_0.size()-1] = 0.0;
  return y3;
}

#endif

