#include <algorithm>
#include <iostream>

template <typename Range>
Range exts( const Range& sexpb, const Range& exb, const Range& beta ){
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

  Range sex ( 2 * sexpb.size() + 1, 0.0 );
  std::reverse_copy(std::begin(sexpb), std::end(sexpb), std::begin(sex) );
  unsigned int k = (beta[0] <= 1.0e-9) ? beta.size() + 1 : beta.size() + 2;
  sex[k-2] = sexpb[0];

  for ( size_t b = 1; b < beta.size(); ++b ){
    sex[k-1] = sexpb[b]*exb[b];
    k += 1;
  }
  return sex;
}
