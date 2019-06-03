#ifndef LEAPR_SIMPLE_CONTINUOUS_INTEGRATE
#define LEAPR_SIMPLE_CONTINUOUS_INTEGRATE


template <typename V>
auto integrate(const V& x, const V& y){
  double value = 0.0;
  value += y[0]*0.5*(x[1]-x[0]);
  for (size_t i = 1; i < x.size()-1; ++i){
    value += y[i]*0.5*(x[i]-x[i-1]);
    value += y[i]*0.5*(x[i+1]-x[i]);
  }
  value += y[x.size()-1]*0.5*(x[x.size()-1]-x[x.size()-2]);
  return value;
}

#endif 
