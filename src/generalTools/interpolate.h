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
  if ( abs(x - x_L) < 1e-8 ){ return y[i]; }
  Float delta = betaGrid[i+1]-betaGrid[i];
  return y[i] + (x-x_L) * (y[i+1]-y[i]) / delta;
}



template <typename Float>
auto terp1( const Float& x1, const Float& y1, const Float& x2, const Float& y2, 
  const Float& x, int i ){
  /* Interpolate one point, where  (x1,y1) and (x2,y2) are the end points of 
   * the line,  (x,y) is the interpolated point, i is the interpolation code,
   * thr6 is the kinematic threshold for i = 6 (thr6 > 0)
   */

  Float a, b, t;

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
  //else if (i == 5) {
  else {
    if (y1 == 0.0) {
      //y = y1;
      return y1;
    } else {
      //y = y1 * exp(log(x/x1) * log(y2/y1) / log(x2/x1));
      return y1 * exp(log(x/x1) * log(y2/y1) / log(x2/x1));
    }
  }

  /* This is commmented out because I have no idea what thr6 is and 
   * I doubt I will anytime soon. Shame
  // coulomb penetrability law (charged particles only)
  else if (i == 6) {
    if (y1 == 0.0) {
      y = y1;
    else
      t = sqrt(x1-thr6);
      b = log((x2*y2)/(x1*y1));
      b = b / (1/t-1/sqrt(x2-thr6));
      a = exp(b/t) * x1 * y1;
      y = (a/x) * exp(-b/sqrt(x-thr6));
    }
  }
  */
}




#endif


