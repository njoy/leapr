#include "generalTools/interpolate.h"
#include "generalTools/print.h"
#include <range/v3/all.hpp>
#include <cmath>


double terps( const std::vector<double>& sd, int nsd, const double& delta, 
  const double& x_val ){
  // Because we need another interpolation function
	
  double yL, yR, y_val, min = -225.0;

  // Index of desired x value
  int xi = x_val/delta;
  // If out of bounds, return 0.0
  if ( xi < 0 or xi > nsd - 1 ){ return 0.0; }
    
  double xL = xi * delta;
  // Don't take the log of a negative number
  yL = sd[xi]   <  0.0 ? min : log( sd[xi]   );
  yR = sd[xi+1] <= 0.0 ? min : log( sd[xi+1] );
    
  y_val = yL + ( x_val - xL ) * ( yR - yL ) / ( delta );

  return y_val >= min ? exp(y_val) : 0.0;
}



template <typename Range, typename Float>
auto terps2( Range sd, int nsd, const Float delta, const Float x_val ){

  std::vector<double> x { 0, 1, 2, 3 }, y { 2, 3, 0, 10 };
  auto xyRange = ranges::view::zip(x,y);
  std::cout << interpolate(xyRange,2.4) << std::endl;

  //auto xVals = ranges::view::linear_distribute(0.0,double(sd.size()*delta),int(sd.size())) | ranges::view::all;
  auto xVals = ranges::view::iota(0,20);//sd.size()*delta,sd.size()) | ranges::view::all;
  auto sd2 = sd | ranges::view::transform([](auto x){return log(x);}) | ranges::view::all;
  //auto xyRange2 = ranges::view::zip(xVals,sd2);
  //std::cout << interpolate(xyRange2,2.4) << std::endl;
  printRange(xVals);
  printRange(sd2);
  std::cout << xVals[0] << std::endl;
  std::cout << sd2[0] << std::endl;
 

  return;
  std::cout << x_val << sd.size() << nsd << delta << std::endl;
  /*
  using std::log;
  //auto xVals = ranges::view::linear_distribute(0.0,sd.size()*delta,sd.size());
  //auto sd2 = sd | ranges::view::transform([](auto x){return log(x);}) | ranges::view::all;
  auto a = ranges::view::iota(0,20);
  auto b = ranges::view::iota(0,20);

  //printRange(xVals);
  //printRange(sd2);
  //return interpolate( ranges::view::zip(xVals,sd2), x_val );
  return interpolate( ranges::view::zip(a,b), x_val );
  std::cout << nsd << delta << x_val << sd.size()<< std::endl;
  */

}



