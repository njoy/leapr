#include <iostream>
#include <vector>

auto bt ( int j, double& pj, const double& x ){
  /* Statistical weight factor
   * for cold hydrogen or deuterium calculation
   */
  int i, k;
  double yy, a, b;
  
  yy = 0.5 * (j+1) * (j+2);
  a = (2*j+3) * exp(-yy*x);
  b = 0;
  for ( auto i = 0; i < 10; ++i ){
    k = 2 * i;
    if (j%2==0) { k = k + 1; }
    yy = 0.5 * k * (k+1);
    b = b + (2*k+1) * exp(-yy*x);
  }
  pj = a / (2 * b);
}

