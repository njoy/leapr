#ifndef SIMPLE_CONTINUOUS_INTERPOLATE
#define SIMPLE_CONTINUOUS_INTERPOLATE


template <typename V, typename F>
F interpolate(const V& xList, const V& yList, const F& x){
  size_t len = xList.size();
  for (size_t i = 0; i < len-1; ++i){
    if (xList[i] <= x and x < xList[i+1]){
      F m = (yList[i+1]-yList[i]) / (xList[i+1]-xList[i]);
      F b = yList[i] - m * xList[i];
      return m*x+b;
    }
  }
  if (xList[xList.size()-1] < x){ return 0.0; }
  return (std::abs(x-xList[len-1])/xList[len-1] < 1e-6) ? yList[len-1] : 0.0;
}

#endif 

