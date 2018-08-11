#include "jprimeLoop_util/sumh.h"
//#include "../../../discre/discre_util/sint.h"

template <typename F, typename A>
//auto jPrime( F& total, int j, const F& be, const F& x, const F& sumConst, 
//  const F& pj, int jj, const A& bex, const A& rdbex, const A& sex, const A& betan, 
//  const F& al, const F& wt, const F& tbart, const F& y, int nbx, bool odd, bool free ){

auto jPrime( F total, int j, F be, F x, F sumConst, 
  F pj, int jj, A bex, A rdbex, A sex, A betan, 
  F al, F wt, F tbart, F y, int nbx, bool odd, bool free ){
  //--sum over the odd or even values of j-prime
  
  // If we want odd values, 1-->9, if even 2-->8.
  auto jpRange = ranges::view::iota(0,5) 
               | ranges::view::transform([odd](auto i){ 
                  return odd ? i*2+1 : i*2; });

  auto bnRange = jpRange | ranges::view::transform([be,j,x](auto jPrime){ 
      return be + ( -j*(j+1) + jPrime*(jPrime+1) ) * x * 0.5; });

  auto tmpRange = jpRange | ranges::view::transform([pj,sumConst,j,y]
      (auto jPrime) { return (2*jPrime+1)*pj*sumConst*4*sumh(j,jPrime,y); });
  // Here, sumConst is equal to 
  //       4*pi/sigma_b * A (if this is a sum over even values) 
  //                             or
  //       4*pi/sigma_b * B (if this is a sum over odd values) 
  // The sumh function is used to calculate the sum over l in Eq. 567-568,
  // which includes the j_l^2(y) term that is the spherical Bessel function
  // of order l

  auto addRange = bnRange | ranges::view::transform([al,wt,free,bex,sex,rdbex,
                                                     betan,tbart,nbx](auto bn){
    if ( free ) {
      // If molecular translations are assumed to be free, we calculate the 
      // S_f(a,b) by using Eq. 569-570. This will be subsequently used in 
      // Eq. 567-568.
      auto ex = -std::pow(al*wt-std::abs(bn),2)/(4*al*wt);
      return (bn > 0.0) ? exp(ex-bn) / sqrt(4*M_PI*al*wt) : 
                          exp(ex)    / sqrt(4*M_PI*al*wt) ;
    }
    else{
      // If the molecular translations are not assumed to be free, we have to
      // calculate S_f(a,b) ourselves, which brings us to use the sint
      // interpolation function we used for discre. 
      return sint(bn,bex,rdbex,sex,betan,betan.size()-1,al,wt,tbart,nbx);
    } } );

  auto summingRange = ranges::accumulate( 
         ranges::view::zip( addRange, tmpRange ),
         std::make_pair(0.0,total),
         [jj](auto l, auto r){
           F add_L = std::get<0>(l); F add_R = std::get<0>(r);
           F tmp_L = std::get<1>(l); F tmp_R = std::get<1>(r);
           return (tmp_R >= 1.0e-6 and jj == 0) ? 
             std::make_pair( add_L + add_R*tmp_R, tmp_L + tmp_R ) :
             std::make_pair( add_L + add_R*tmp_R, tmp_L); 
           }
      );
  total = std::get<1>(summingRange);
  return  std::get<0>(summingRange);
}
