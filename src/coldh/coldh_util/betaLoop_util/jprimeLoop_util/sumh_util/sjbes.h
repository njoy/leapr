#include <iostream>
#include <vector>


double sjbes( int n, double y ){
  /* Calculates spherical bessel functions for cold hydrogen or deuterium 
   * calculation. These will be used in Eq. 567 - 568
   * The spherical bessel function is of nth order, evaluated at position y.
   * In the case of cold hydrogen and deuterium calculations, y is defined by
   */

  int k, l, nm;
  double w, sj, t1, t2, t3;

  // check for large or negative arguments
  if ( n >= 3.0e4 or y > 3.0e4 or y < 0.0 or n < 0.0 ){
    return 0.0;
  }

  // compute normal values
  if (y <= 7.0e-4) {
    if      (n  == 0) { return 1; } 
    else if (n  > 10) { return 0; } 
    else {
      t3 = 1;
      for ( int i = 3; i < 2*n + 3; i = i + 2 ){
        t3 = t3 * i;
      }
      return std::pow(y,n)/t3;
    }
  }
  else {
    w = y < 0.2 ? 1 - y*y * ( 1 - y*y/20 ) / 6 : 
                  sin(y) / y;
    if (n == 0) {
      return w;
    } 
    else {
      if      (y >= 100.0) { l = int( y/50 + 18 ); } 
      else if (y >= 10.0 ) { l = int( y/10 + 10 ); } 
      else if (y >  1.0  ) { l = int( y/2  + 5  ); } 
      else                 { l = 5; }

      nm = y > n ? y + l : n + l;
      t3 = 0;
      t2 = 2.0e-38;
      for ( auto i = 0; i < nm; ++i ){
        k = nm - i - 1;
        t1 = (2*k + 3) * t2 / y - t3;

        if (n == k){ sj = t1;}
        if (std::abs(t1) >= 1.0e25) {
          t1 = t1 / 1.0e25;
          t2 = t2 / 1.0e25;
          sj = sj / 1.0e25;
        }
        t3 = t2;
        t2 = t1;
      }
      return w * sj / t1;
    }
  }

}
