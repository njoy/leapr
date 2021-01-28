#include <iostream>
#include <vector>
#include "jprimeLoop_util/sumh.h"
#include <cmath>

template <typename Float, typename Range>
auto jPrime( int j, const Float& be, const Float& x, const Float& sumConst, 
  const Float& pj, const Range& bex, const Range& rdbex, const Range& sex, 
  const Range& betan, const Float& alphaWgt, const Float& tbart, const Float& y, 
  int nbx, bool odd, bool free){
  //--sum over the odd or even values of j-prime
  
  Float add, snl = 0, tmp, bn, ex;
  int start, end;

  // Select start and end values for loop. If we want odd values, 1-->9,
  // if even 2-->8.
  if ( odd == true ){ start = 1; end = 10; }
  else              { start = 0; end = 9;  }

  for ( auto jPrime = start; jPrime <= end; jPrime = jPrime + 2 ){

    bn = be + (-j*(j+1) + jPrime*(jPrime+1)) * x * 0.5;

    tmp = (2*jPrime+1) * pj * sumConst * 4.0 * sumh(j,jPrime,y);
    // Here, sumConst is equal to 
    //       4*pi/sigma_b * A (if this is a sum over even values) 
    //                             or
    //       4*pi/sigma_b * B (if this is a sum over odd values) 
    // The sumh function is used to calculate the sum over l in Eq. 567-568,
    // which includes the j_l^2(y) term that is the spherical Bessel function
    // of order l

    if ( free ) {
      // If molecular translations are assumed to be free, we calculate the 
      // S_f(a,b) by using Eq. 569-570. This'll be used in Eq. 567-568.
      ex = -std::pow(alphaWgt-std::fabs(bn),2)/(4*alphaWgt);
      if ( bn > 0.0 ){ ex -= bn; }
      add = exp(ex)/sqrt(4*M_PI*alphaWgt);
    }
    else{
      // If the molecular translations are not assumed to be free, we have to
      // calculate S_f(a,b) ourselves, which brings us to use the sint
      // interpolation function we used for discre. 
      add = sint(bn,bex,rdbex,sex,betan,betan.size()-1,alphaWgt,tbart,nbx);
    }
    snl += tmp * add;
  }
  return snl;

}
