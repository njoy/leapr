#include <range/v3/all.hpp>
#include <iostream>
#include "generalTools/print.h"
#include <algorithm>

template <typename Range, typename Float>
auto search( Range xRange, Float x, int i, int left, int right ){
  if ( xRange[i] <= x and x <= xRange[i+1] ){ return i; }
  if ( x > xRange[i] ){ left  = i; }
  else                { right = i; }
  i = (left + right)*0.5;
  return search(xRange,x,i,left,right);
}


template <typename RangeZip, typename Float>
Float interpolate( RangeZip xyRange, Float x ){
  int len = xyRange.size();
  auto xVec = xyRange | ranges::view::keys;
  auto yVec = xyRange | ranges::view::values;
  if ( x < std::get<0>(xyRange[0]) or 
       x > std::get<0>(xyRange[len-1]) ){
    return 0.0;
  }

  int index = search(xVec, x, int(len*0.5), 0, len);
  Float b = yVec[index];
  Float m = (yVec[index+1]-yVec[index])/(xVec[index+1]-xVec[index]);
  return m*(x-xVec[index])+b;
}





