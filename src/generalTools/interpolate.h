#ifndef LEAPR_GENERAL_TOOLS_INTERPOLATE
#define LEAPR_GENERAL_TOOLS_INTERPOLATE

#include <range/v3/all.hpp>
#include <iostream>
#include <algorithm>
#include <tuple>

template <typename Range, typename Float>
auto search( Range xRange, Float x, int i, int left, int right ){
  if ( xRange[i] <= x and x <= xRange[i+1] ){ return i; }
  if ( x > xRange[i] ){ left  = i; }
  else                { right = i; }
  i = (left + right)*0.5;
  return search(xRange,x,i,left,right);
}


template <typename RangeZip, typename Float>
Float interpolate( RangeZip&& xyRange, Float& x, Float&& defaultVal=0.0 ){
  int len = xyRange.size();
  auto xVec = xyRange | ranges::view::keys;
  auto yVec = xyRange | ranges::view::values;
  if ( x < xVec[0] or x > xVec[len-1] ){ return defaultVal; }
  int index = search(xVec, x, int(len*0.5), 0, len);
  Float b = yVec[index];
  Float m = (yVec[index+1]-yVec[index])/(xVec[index+1]-xVec[index]);
  return m*(x-xVec[index])+b;
}

template <typename RangeZip, typename Float>
Float interpolateLog( RangeZip&& xyRange, Float& x ){
  using std::log; using std::exp;
  auto xVec = xyRange | ranges::view::keys;
  auto yVec = xyRange | ranges::view::values
                      | ranges::view::transform([](auto x){return log(x);});
  return exp(interpolate( ranges::view::zip(xVec,yVec), x ));
}


template <typename Float, typename Range>
Float interpolate( const Range& y, const Float& x, const Range& betaGrid ){
  if ( x >= betaGrid[betaGrid.size()-1] ){ return 0.0; }
  if ( x < betaGrid[0] ){ return 0.0; }
 
  unsigned int i = 0;
  Float x_L = 0.0;
  for ( size_t j = 0; j < betaGrid.size(); ++j ){
    if ( x < betaGrid[j+1] ){ 
      i = j;
      x_L = betaGrid[j];
      break;
    }
  }
  if ( abs(x - x_L) < 1e-8 ){ return y[i]; }
  Float delta = betaGrid[i+1]-betaGrid[i];
  return y[i] + (x-x_L) * (y[i+1]-y[i]) / delta;
}





#endif


