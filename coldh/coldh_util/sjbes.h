#include <iostream>
#include <vector>


double sjbes( int n, double x ){
  /* Bessel functions for cold hydrogen or deuterium calculation
   */
  int i, k, l, iii, kmax, nm;
  double w, bessel, y, z, sj, t1, t2, t3, sjbes;

  // check for large or negative arguments
  if ( n >= 3.0e4 or x > 3.0e4 or x < 0.0 or n < 0.0 ){
    return 0.0;
  }

  // compute normal values
  if (x <= 7.0e-4) {
    if      (n  == 0) { return 1; } 
    else if (n  > 10) { return 0; } 
    else {
      t3 = 1;
      for ( int i = 3; i < 2*n + 3; i = i + 2 ){
        t3 = t3 * i;
      }
      return std::pow(x,n)/t3;
    }
  }
  else {
    w = x < 0.2 ? 1 - x*x * ( 1 - x*x/20 ) / 6 : 
                  sin(x) / x;
    if (n == 0) {
      return w;
    } 
    else {
      if      (x >= 100.0) { l = int( x/50 + 18 ); } 
      else if (x >= 10.0 ) { l = int( x/10 + 10 ); } 
      else if (x >  1.0  ) { l = int( x/2  + 5  ); } 
      else                 { l = 5; }

      nm = x > n ? x + l : n + l;
      t3 = 0;
      t2 = 2.0e-38;
      for ( auto i = 0; i < nm; ++i ){
        k = nm - i - 1;
        t1 = (2*k + 3) * t2 / x - t3;

        if (n == k){ sj = t1;}
        if (abs(t1) >= 1.0e25) {
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
