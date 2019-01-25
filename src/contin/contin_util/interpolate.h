#include <iostream>

template <typename F, typename A>
double interpolate( const A& y, const F& delta2, const F& x, const A& betaGrid ){
  /* Inputs
   * ------------------------------------------------------------------------
   * y     : evenly spaced vector to be interpolated
   * delta : spacing in y
   * x     : desired value
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Interpolate
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * interpolated value
   */

  // This is the index just to the left of desired x point
  unsigned int i = x / delta2; 
  if ( x >= (y.size()-1)*delta2 ){ return 0.0; }
  double x_L = i * delta2; 
                     
  // Return 0.0 if out of range (above), else return interpolated value
  return ( i+2 > y.size() ) ? 0.0 : y[i] + (x-x_L) * (y[i+1]-y[i]) / delta2;
  std::cout << betaGrid.size() << std::endl;

  /*
  if ( x >= betaGrid[betaGrid.size()-1] ){ return 0.0; }
  if ( x < betaGrid[0] ){ return 0.0; }
 
  // This is the index just to the left of desired x point
  unsigned int i = 0;
  F x_L = 0.0;
  for ( size_t j = 0; j < betaGrid.size(); ++j ){
    if ( x < betaGrid[j+1] ){ 
      i = j;
      x_L = betaGrid[j];
      break;
    }
  }
  if ( abs(x - x_L) < 1e-8 ){ return y[i]; }

  F delta = betaGrid[i+1]-betaGrid[i];
                     
  // Return 0.0 if out of range (above), else return interpolated value
  return y[i] + (x-x_L) * (y[i+1]-y[i]) / delta;
  std::cout << delta2 << std::endl;
  */


}


