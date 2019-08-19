#include "generalTools/interpolate.h"
#include "generalTools/print.h"
#include <range/v3/all.hpp>

template <typename Range, typename Float>
auto terps( Range sd, const Float delta, const Float x_val, int nsd ){
  if ( x_val < 0.0 or x_val > nsd*delta ){ return 0.0; }
  auto xVals = ranges::view::iota(0,int(sd.size())) 
             | ranges::view::transform([delta](auto x){return delta*x;});
  return interpolateLog(ranges::view::zip(xVals,sd),x_val);
}



