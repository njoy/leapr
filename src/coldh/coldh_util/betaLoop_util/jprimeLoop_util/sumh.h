#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/cn.h"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/sjbes.h"
#include <range/v3/all.hpp>

auto sumh(int j, int jp, double y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation, which follows Eq. 567 - 568
   */ 
  if      (j  == 0) { return std::pow( sjbes(jp, y) * cn(j, jp, jp), 2 ); } 
  else if (jp == 0) { return std::pow( sjbes(j , y) * cn(j, 0,  j ), 2 ); }

  return ranges::accumulate(
           ranges::view::iota(std::abs(j-jp)+1,std::min(j+jp+2,std::abs(j-jp)+11)) 
         | ranges::view::transform( [y,j,jp](auto n){ 
             return std::pow( sjbes(n-1,y) * cn(j,jp,n-1), 2 ); }),
         0.0);

}
