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


  // Here we reverse the vector sexpb and put it in the beginning of sex.
  
  unsigned int k = beta[0] <= 1.0e-9 ? beta.size() + 1 : beta.size() + 2;

  auto zeros1 = ranges::view::iota( 1, int(beta.size()+2)) 
              | ranges::view::transform([](auto){return 0.0;});

  auto zeros2 = ranges::view::iota( 1, int(beta.size()-k+4) )
              | ranges::view::transform( [](auto){ return 0.0; } );

  // { sexpb[1]*exb[1]   sexpb[2]*exb[2]   sexpb[3]*exb[3] .... }
  auto sexpb_times_exb = ranges::view::zip(sexpb,exb) 
                       | ranges::view::transform([](auto t){
                           return std::get<0>(t)*std::get<1>(t); }) 
                       | ranges::view::slice(1,ranges::end);

  // First k-2 entries of { sexpb[n], ... , sexpb[0], 0, 0, ... 0 }
  auto sex = ranges::view::concat( sexpb | ranges::view::reverse, zeros1 )
           | ranges::view::slice( 0, int(k-2) );
 
  // Potentially making two s1's (in the example at the beginning of this 
  // function)
  return ranges::view::concat( 
           sex,
           ranges::view::single(sexpb[0]),
           sexpb_times_exb,
           zeros2);
  /*                |--sex--|   single |-- sexpb*exb --|  zero(s)
   *              [ s3, s2, s1,   s1,  s2*e2*e2, s3*e3*e3,   0   ]
   *                              or
   *                |--sex--|   single |-- sexpb*exb --|  zero(s)
   *              [ s3, s2, s1,        s2*e2*e2, s3*e3*e3, 0, 0 ]
   */
 

  // sex --> ssm_i * exp( -beta ) 
  // S(a,b) = exp( -beta ) * S(a,-b)    Eq. 509
  // Notice that we only apply this to half of the sex vector, since only
  // half of it needs to be flipped
}
