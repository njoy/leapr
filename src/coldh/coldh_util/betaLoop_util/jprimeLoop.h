#include <iostream>
#include <vector>
#include "jprimeLoop_util/sumh.h"
//#include "../../../discre/discre_util/sint.h"

template <typename F, typename I, typename A>
auto jPrime( F& total, I j, const F& be, const F& x, const F& sumConst, 
  const F& pj, I jj, const A& bex, const A& rdbex, const A& sex, const A& betan, 
  const F& al, const F& wt, const F& tbart, const F& y, I nbx, bool odd, bool free ){

  //--sum over the odd or even values of j-prime
  
  F add, snl = 0, pi = 3.14159265358979, tmp, bn, ex;
  I start, end;

  // Select start and end values for loop. If we want odd values, 1-->9,
  // if even 2-->8.
  if ( odd == true ){ start = 1; end = 10; }
  else              { start = 0; end = 9;  }
  auto jPrimeRange = ranges::view::iota(0,5) 
              | ranges::view::transform([odd](auto i){ 
                  return odd ? i*2 + 1 : i*2; });
  auto bnRange = jPrimeRange | ranges::view::transform([be,j,x](auto jPrime){ 
      return be + ( -j*(j+1) + jPrime*(jPrime+1) ) * x * 0.5; });

  auto tmpRange = jPrimeRange | ranges::view::transform([pj,sumConst,j,y]
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

  auto snlRange = 
    ranges::accumulate( 
      ranges::view::zip( addRange, tmpRange ) 
    | ranges::view::transform( [](auto t){ 
        return std::make_pair(std::get<0>(t)*std::get<1>(t),std::get<1>(t));}), 
    std::make_pair(0.0,total),
    [jj](auto l, auto r){
    F snl_L = std::get<0>(l); F snl_R = std::get<0>(r);
    F tmp_L = std::get<1>(l); F tmp_R = std::get<1>(r);
    if ( tmp_R >= 1.0e-6 and jj == 0 ){ 
      return std::make_pair(snl_L+snl_R,tmp_L+tmp_R); 
    }
    return std::make_pair(snl_L+snl_R,tmp_L); 
    });
  //std::cout << std::get<0>(snlRange) << "     " << std::get<1>(snlRange) << std::endl;
  total = std::get<1>(snlRange);
  return std::get<0>(snlRange);

  //total = std::get<1>(snlRange);
  //return std::get<0>(snlRange);
            
  /*
  auto test1 = ranges::view::iota(1,5);
  auto test2 = ranges::view::iota(5,9);
  auto test3 = ranges::view::zip(test1,test2);
  std::cout << test1 << std::endl;
  std::cout << test2 << std::endl;
  auto sum1 = ranges::accumulate(test3,std::make_pair(0.0,0.0),[](auto l,auto r){
    return std::make_pair(std::get<0>(l)+std::get<0>(r),std::get<1>(l)+std::get<1>(r));
  });
  std::cout << std::get<0>(sum1) << "    " << std::get<1>(sum1) << std::endl;
  */

//  std::cout << snlRange << std::endl;
/*


  for ( auto jPrime = start; jPrime < end; jPrime = jPrime + 2 ){

    bn = be + ( -j*(j+1) + jPrime*(jPrime+1) ) * x * 0.5;

    tmp = (2*jPrime+1) * pj * sumConst * 4 * sumh(j,jPrime,y);
    // Here, sumConst is equal to 
    //       4*pi/sigma_b * A (if this is a sum over even values) 
    //                             or
    //       4*pi/sigma_b * B (if this is a sum over odd values) 
    // The sumh function is used to calculate the sum over l in Eq. 567-568,
    // which includes the j_l^2(y) term that is the spherical Bessel function
    // of order l

    if (jj == 0 and tmp >= 1.0e-6) { total += tmp; }

    if ( free ) {
      // If molecular translations are assumed to be free, we calculate the 
      // S_f(a,b) by using Eq. 569-570. This will be subsequently used in 
      // Eq. 567-568.
      ex = -std::pow(al*wt-std::abs(bn),2)/(4*al*wt);
      add = (bn > 0.0) ? exp(ex-bn) / sqrt(4*pi*al*wt) : 
                         exp(ex)    / sqrt(4*pi*al*wt) ;
    }
    else{
      // If the molecular translations are not assumed to be free, we have to
      // calculate S_f(a,b) ourselves, which brings us to use the sint
      // interpolation function we used for discre. 
      add = sint(bn,bex,rdbex,sex,betan,betan.size()-1,al,wt,tbart,nbx);
    }

    snl += tmp * add;

  }
  std::cout << snl << "     " << total << std::endl;
  return snl;

  */
  /*
  */
}
