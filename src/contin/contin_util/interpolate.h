#include <iostream>

template <typename Float, typename Range>
Float interpolate( const Range& y, const Float& x, const Range& betaGrid ){
  /* Inputs
   * ------------------------------------------------------------------------
   * y     : vector to be interpolated
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

  if ( x >= betaGrid[betaGrid.size()-1] ){ return 0.0; }
  if ( x < betaGrid[0] ){ return 0.0; }
 
  // This is the index just to the left of desired x point
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


