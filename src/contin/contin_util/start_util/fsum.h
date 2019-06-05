#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH
#include <range/v3/all.hpp>
#include <tuple>
#include <iostream>
#include "generalTools/trapezoidIntegral.h"

using std::cosh; using std::pow;

template <typename Float, typename Tuple>
auto even(Float tau, int n, Tuple xyPair){
  Float beta = std::get<0>(xyPair);
  Float pVal = std::get<1>(xyPair);
  return 2.0*pVal*cosh(beta*tau)*pow(beta,n);
}

template <typename Float, typename Tuple>
auto odd(Float tau, int n, Tuple xyPair){
  Float beta = std::get<0>(xyPair);
  Float pVal = std::get<1>(xyPair);
  return 2.0*pVal*sinh(beta*tau)*pow(beta,n);
}

template <typename Range, typename Float>// = Range::value_type >
auto fsum( int n, Range beta_P_zipped, Float tau ){
  /* Inputs
   * ------------------------------------------------------------------------
   * n       : appears in equation being evaluated
   * p       : P(beta), which is defined in Eq. 507
   * tau     : appears in equation being evaluated
   * betas   : beta grid for integral
   *
   * Operations
   * ------------------------------------------------------------------------
   * * Compute an integral of the form
   *   int 0->infty  2 * P(beta) * beta^n * sinh( tau*beta ) dbeta
   *                                or
   *   int 0->infty  2 * P(beta) * beta^n * cosh( tau*beta ) dbeta
   *  depending on whether n is odd of even, respectively.
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * The computed integral is returned
   *
   */

  if (n%2 == 0){ 
    auto evenLambda = [tau, n](auto xyPair){ return even(tau,n,xyPair); };
    return trapezoidIntegral( beta_P_zipped, evenLambda );
  }
  else { 
    auto oddLambda  = [tau, n](auto xyPair){ return odd(tau,n,xyPair);  };
    return trapezoidIntegral( beta_P_zipped, oddLambda );
  }

 
}







#endif
