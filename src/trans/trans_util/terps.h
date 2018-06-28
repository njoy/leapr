

template <typename floatT, typename arrayT>
floatT terps( const arrayT& sd, int nsd, const floatT& delta, 
  const floatT& xVal ){
  // Because we need another interpolation function
	
  floatT yL, yR, yVal, min = -225.0;

  // Index of desired x value
  int xi = xVal/delta;

  // If out of bounds, return 0.0
  if ( xi < 0 or xi > nsd - 1 ){ return 0.0; }
    
  floatT xL = xi * delta;

  // Don't take the log of a negative number
  yL = sd[xi]   <  0.0 ? min : log( sd[xi]   );
  yR = sd[xi+1] <= 0.0 ? min : log( sd[xi+1] );
    
  yVal = yL + ( xVal - xL ) * ( yR - yL ) / ( delta );

  return yVal >= min ? exp(yVal) : 0.0;
}
