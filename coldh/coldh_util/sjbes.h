#include <iostream>
#include <vector>


auto sjbes( int n, double x ){
  /* Bessel functions for cold hydrogen or deuterium calculation
   */
  int i, k, l, iii, kmax, nm;
  double w, bessel, y, z, sj, t1, t2, t3, sjbes;

  // check for large arguments
  if (n >= 30000 or x > 3.0e4) {
      sjbes = 0;
      return sjbes;
  }

  // check for bad arguments
  if (x < 0.0 or n < 0) {
    sjbes = 0;
    return sjbes;
  }

  // compute normal values
  if (x <= 7.0e-4) {
    w = 1;
    if (n == 0) {
      bessel = w;            
    } else if (n > 10) {
      bessel = 0;
    } else {
      t1 = 3;
      t2 = 1;
      for ( int i = 0; i < n; ++i ){
        t3 = t2 * x / t1;
        t1 = t1 + 2;
        t2 = t3;
      }
      bessel = t3;
    }
  } else {
    if (x < 0.2) {
      y = x * x;
      w = 1 - y * ( 1 - y / 20 ) / 6;
    } else {
      w = sin(x) / x;
    }
    if (n == 0) {
      bessel = w;
    } else {
      if (x >= 100.0) {
        l = int(x/50+18);
      } else if (x >= 10.0) {
        l = int(x/10+10);
      } else if (x > 1.0) {
        l = int(x/2+5);
      } else {
        l = 5;
      }
      iii = int(x);
      kmax = n;
      if (iii > n) kmax = iii;
      nm = kmax + l;
      z = 1 / x;
      t3 = 0;
      t2 = 2.0e-38;
      for ( auto i = 0; i < nm; ++i ){
        k = nm - i - 1;
        t1 = (2*k + 3) * z * t2 - t3;
        if (n == k){ sj = t1;}
        if (abs(t1) >= 1.0e25) {
          t1 = t1 / 1.0e25;
          t2 = t2 / 1.0e25;
          sj = sj / 1.0e25;
          
      std::cout << "HERE"<< std::endl;
        }
        t3 = t2;
        t2 = t1;
      }
      bessel = w * sj / t1;
    }
  }
  sjbes = bessel;
  return sjbes;

}
