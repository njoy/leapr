#include "simple/continuous/interpolate.h"
#include <iostream>


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

