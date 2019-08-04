#include "coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/cn.h"
#include <boost/math/special_functions/bessel.hpp>
#include <range/v3/all.hpp>

template <typename Float>
auto sumh(int j, int jp, Float y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation, which follows Eq. 567 - 568
   */ 
  using boost::math::sph_bessel;
  if      (j  == 0) { return std::pow( sph_bessel(jp, y) * cn(j, jp, jp), 2 ); } 
  else if (jp == 0) { return std::pow( sph_bessel(j , y) * cn(j, 0,  j ), 2 ); }

  return ranges::accumulate(
           ranges::view::iota(std::abs(j-jp)+1,std::min(j+jp+2,std::abs(j-jp)+11)) 
         | ranges::view::transform( [y,j,jp](auto n){ 
             return std::pow( sph_bessel(n-1,y) * cn(j,jp,n-1), 2 ); }),
         0.0);

}
