
#ifndef THERMR_LEGNDR_HH
#define THERMR_LEGNDR_HH
inline auto legndr( double x, std::vector<double>& p, int np ){
 /*--------------------------------------------------------------------
  * Generate Legendre polynomials at x by recursion.
  * Place p(subl) in p(l+1).
  *--------------------------------------------------------------------
  */

  // Defined to be the first two Legendre polynomial values.
  p[0] = 1; p[1] = x;

  if (np > 1){
    // If you need the 2nd or higher order, must actually calculate things
    for ( size_t i = 0; int(i) < np - 1 and i < p.size()-2; ++i ){
      p[i+2] = 2*x*p[i+1] - p[i] - (x*p[i+1]-p[i]) / (i+2);
    }
  }
}

#endif
