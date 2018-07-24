#include <algorithm>
#include <iostream>
#include <range/v3/all.hpp>

template <typename A>
A exts( const A& sexpb, const A& exb, const A& beta ){
  /* When used by discre, sexpb is a vector populated with sym_sab entries for 
   * a given alpha and temp and increasing values of beta. Also exts is 
   * a vector populated with exp( -beta * sc / 2 ) entries. The purpose for 
   * the exts vector is to change S(a,b) --> S(a,-b) by multiplying it by 
   * exp( -beta ). 
   *
   *          Given beta = [ b1, b2, b3 ], sexpb = [ s1, s2, s3 ], 
   *          and exb = [ e1, e2, e3, e4, e5, e6, e7 ], then sex is
   *
   *              [ s3, s2, s1, s1, s2*e2*e2, s3*e3*e3, 0 ]
   *                          if b1 > 1e-9,
   *                       
   *              [ s3, s2, s1, s2*e2*e2, s3*e3*e3, 0, 0 ]
   *                               else
   */

  A sex ( 2 * sexpb.size() + 1, 0.0 );

  // Here we reverse the vector sexpb and put it in the beginning of sex.
  
  unsigned int k = beta[0] <= 1.0e-9 ? beta.size() + 1 : beta.size() + 2;

  std::reverse_copy(std::begin(sexpb), std::end(sexpb), std::begin(sex) );
  auto sexRange = ranges::view::concat(
                    sexpb | ranges::view::reverse, 
                    ranges::view::iota(1,int(2*sexpb.size()+1-sexpb.size()+1)) 
                  | ranges::view::transform([](auto){return 0.0;}) );

  auto sexRange2 = ranges::view::concat(
                     sexRange|ranges::view::slice(0,int(k-2)),
                     ranges::view::single(sexpb[0]));
                     //sexRange|ranges::view::slice(int(k-1),ranges::end));

  sex[k-2] = sexpb[0];


  auto sexpb_exb = ranges::view::zip(sexpb,exb) | ranges::view::transform([](auto t){return std::get<0>(t)*std::get<1>(t); }) | ranges::view::slice(1,ranges::end);

  auto finalSexRange = ranges::view::concat(sexRange2,sexpb_exb);
  auto finalSexRange2 = ranges::view::concat(finalSexRange,ranges::view::iota(1,int(2*sexpb.size()+1-finalSexRange.size()+1))|ranges::view::transform([](auto){return 0.0;}));

  for ( size_t b = 1; b < beta.size(); ++b ){
    std::cout << sexpb[b]*exb[b] << std::endl;
    // sex --> ssm_i * exp( -beta ) 
    // S(a,b) = exp( -beta ) * S(a,-b)    Eq. 509
    // Notice that we only apply this to half of the sex vector, since only
    // half of it needs to be flipped
    sex[k-1] = sexpb[b]*exb[b];
    k += 1;
  }
  std::cout << (sex|ranges::view::all) << std::endl;
  std::cout << std::endl;
  //std::vector<double> outVec = finalSexRange2|ranges::to_vector;
  //std::cout << finalSexRange2|ranges::to_vector << std::endl;
  std::cout << std::endl;
  return finalSexRange2;
  //return std::make_tuple(sex,sex);
}
