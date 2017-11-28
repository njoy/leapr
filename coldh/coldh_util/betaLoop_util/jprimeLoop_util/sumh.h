#include <iostream>
#include <vector>
#include "sumh_util/cn.h"
#include "sumh_util/sjbes.h"


auto sumh(int j, int jp, double y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation.
   */ 
  if      (j  == 0) { return std::pow( sjbes(jp, y) * cn(j, jp, jp), 2 ); } 
  else if (jp == 0) { return std::pow( sjbes(j , y) * cn(j, 0,  j ), 2 ); }

  int ipk, n, imk;
  double sum1 = 0;

  imk = abs(j-jp) + 1;
  ipk = std::min(j+jp+2, abs(j-jp)+11);
  for ( auto n = imk; n < ipk; ++n ){
    sum1 += std::pow( sjbes(n-1,y) * cn(j,jp,n-1), 2 );
  }
  return sum1;
}
