#include <vector>

double interpolate( const std::vector<double>& y_vec, 
                    const double& delta, const double& x ){
  /* This takes in some vector y (which has a corresponding evenly spaced
  * x vector, with spacing delta). It finds and returns the linearly 
  * interpolated value approximation at point x
  */

  int i = x / delta; // This is the index just to the left
                     // of my desired x point

  // Return 0.0 if out of range (above) 
  if ( i >= y_vec.size() - 1 ){ return 0.0; }

  // Check to make sure that the desired x point is within range
  double x_L = i * delta; // tabulated x values to the left and 

  return y_vec[i] + ( x - x_L ) * ( y_vec[i+1] - y_vec[i] ) / delta;
}



