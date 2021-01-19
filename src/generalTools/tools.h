#ifndef LEAPR_GENERAL_TOOLS
#define LEAPR_GENERAL_TOOLS


#include <range/v3/all.hpp>
#include <iostream>
#include <algorithm>
#include <tuple>


template <typename Range, typename Callable >
auto trapezoidIntegral( Range inputXY, Callable callable ){
  auto xVec = inputXY | ranges::view::keys;
  //auto yVec = inputXY | ranges::view::values;
  auto binWidths = xVec | ranges::view::sliding(2) 
                        | ranges::view::transform([](auto pair){ 
                            return pair[1]-pair[0]; } );
  auto argument = inputXY | ranges::view::transform(callable);
  auto outputWindows = argument | ranges::view::sliding(2);
  auto trapezoid = [](auto binWidth, auto leftRightPair){ 
    return (leftRightPair[0]+leftRightPair[1])*0.5*binWidth;
  };
  auto integral = ranges::view::zip_with(trapezoid,binWidths,outputWindows);
  return ranges::accumulate(integral,0.0);
}





template <typename Float>
void swap( Float& a, Float& b ){
    Float c = a;
    a = b;
    b = c;
}




template <typename Range, typename Float>
auto search( const Range& xRange, const Float& x, int i, int left, int right ){
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
Float interpolate( const Range& y, const Float& x, const Range& betaGrid, Float defaultVal = 0.0 ){
  if ( x >= betaGrid[betaGrid.size()-1] ){ return defaultVal; }
  if ( x < betaGrid[0] ){ return defaultVal; }
 
  unsigned int i = 0;
  Float x_L = 0.0;
  for ( size_t j = 0; j < betaGrid.size(); ++j ){
    if ( x < betaGrid[j+1] ){ 
      i = j;
      x_L = betaGrid[j];
      break;
    }
  }

  Float delta = betaGrid[i+1]-betaGrid[i];
  return y[i] + (x-x_L) * (y[i+1]-y[i]) / delta;
}


template <typename Float, typename Range>
Float interpolate( const Range& y, const Float& x, const Float& delta, int np,
                   Float defaultVal = 0.0 ){
  if ( x >= (np-1)*delta or x < 0 ){ return defaultVal; }
  unsigned int i = 0;
  auto xGrid = ranges::view::iota(0,int(y.size())) 
             | ranges::view::transform([delta](int i){return delta*i;});
  i = search(xGrid, x, int(y.size()*0.5), 0, y.size());
  return y[i] + (x-delta*i) * (y[i+1]-y[i]) / delta;
}



template <typename Float>
auto terp1( const Float& x1, const Float& y1, const Float& x2, const Float& y2, 
  const Float& x, int i ){
  /* Interpolate one point, where  (x1,y1) and (x2,y2) are the end points of 
   * the line,  (x,y) is the interpolated point, i is the interpolation code,
   * thr6 is the kinematic threshold for i = 6 (thr6 > 0)
   */

  // make sure x2 .ne. x1
  if (x2 == x1) {
    return y1;
  }

  // y is constant
  if (i == 1 or y2 == y1 or x == x1) {
    return y1;
  } 

  // y is linear in x
  else if (i == 2) {
    return y1 + (x-x1) * (y2-y1) / (x2-x1);
  }

  // y is linear in ln(x)
  else if (i == 3) {
    return y1 + log(x/x1) * (y2-y1) / log(x2/x1);
  }

  // ln(y) is linear in x
  else if (i == 4) {
    return y1 * exp((x-x1) * log(y2/y1) / (x2-x1));
  }
  // ln(y) is linear in ln(x)
  else {
    if (y1 == 0.0) {
      return y1;
    } else {
      return y1 * exp(log(x/x1) * log(y2/y1) / log(x2/x1));
    }
  }
}

#endif


