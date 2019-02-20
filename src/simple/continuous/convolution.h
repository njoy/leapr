#include "simple/continuous/interpolate.h"


template <typename V>
V convolve(V x, V y1, V y2){ 
  V y3 (x.size());
  for (size_t i = 0; i < x.size(); ++i ){
    for (size_t j = 0; j < x.size()-1; ++j ){
      auto val_L = y1[j]*interpolate(x,y2,x[i]-x[j]);
      auto val_R = y1[j+1]*interpolate(x,y2,x[i]-x[j+1]);
      y3[i] = (val_L + val_R) * 0.5 * (x[j+1]-x[j]);
    }
  }
  return y3;
}


