#include <iostream>
#include <vector>
#include "sumh_util/cn.h"
#include "sumh_util/sjbes.h"


auto sumh(int j, int jp, double y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation, which follows Eq. 567 - 568
   */ 
  if      (j  == 0) { return std::pow( sjbes(jp, y) * cn(j, jp, jp), 2 ); } 
  else if (jp == 0) { return std::pow( sjbes(j , y) * cn(j, 0,  j ), 2 ); }

  int end, n, start;
  double sum1 = 0;

  start = abs(j-jp) + 1;
  end = std::min(j+jp+2, abs(j-jp)+11);
  for ( auto n = start; n < end; ++n ){
    sum1 += std::pow( sjbes(n-1,y) * cn(j,jp,n-1), 2 );
    // sjbes is used to calculate the spherical bessel function j_l
  }
  return sum1;
}
