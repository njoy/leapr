#include <iostream>
#include <vector>


std::vector<double> exts( const std::vector<double>& sexpb,
  const std::vector<double>& exb, const std::vector<double>& beta ){
  /* When used by discre, sexpb is a vector populated with sym_sab entries for 
   * a given alpha and temp and increasing values of beta.
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

  std::vector<double> sex ( 2 * sexpb.size() + 1, 0.0 );

  // Here we reverse the vector sexpb and put it in the beginning of sex.
  
  std::reverse_copy(std::begin(sexpb), std::end(sexpb), std::begin(sex) );

  int k = beta[0] <= 1.0e-9 ? beta.size() + 1 : beta.size() + 2;
  sex[k-2] = sexpb[0];

  for ( int b = 1; b < beta.size(); ++b ){
    sex[k-1] = sexpb[b]*exb[b]*exb[b];
    k += 1;
  }
  return sex;
}
