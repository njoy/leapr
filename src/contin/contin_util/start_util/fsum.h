
#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH

#include <range/v3/all.hpp>
#include <iostream>

template<typename T, typename V>
T fsum( const int n, const V& p, const T& tau, const T& delta_b ){
  /* Inputs
   * ------------------------------------------------------------------------
   * n       : appears in equation being evaluated
   * p       : P(beta), which is defined in Eq. 507 
   * tau     : appears in equation being evaluated
   * delta_b : increment size used for integral. Doing a Riemann sum with 
   *           rectangle sum delta_b wide
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Compute an integral of the form 
   *
   *      infty
   *   int       2 * P(beta) * beta^n * sinh( tau*beta ) dbeta
   *      0
   *                                or
   *      infty
   *   int       2 * P(beta) * beta^n * cosh( tau*beta ) dbeta
   *      0
   *
   *  depending on whether n is odd or even, respectively.
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * The computed integral is returned
   *
   */

  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to 
                                    // help differentiate betwene sinh and 
                                    // cosh while evaluating the integrand

  // Create a range of beta values, [ 0, delta_b, 2*delta_b, 3*delta_b, ...]
  auto b = ranges::view::zip(
             ranges::view::iota(0,int(p.size()))
           | ranges::view::transform([delta_b](auto x){return delta_b*x;}),
	   p );

  auto funcVal = 
    b | ranges::view::transform([tau,n,even,p,delta_b](auto x){ 
          T p_i = std::get<1>(x), beta = std::get<0>(x);
          T val = even ? 2.0 * p_i * cosh(beta*tau) * std::pow(beta,n) :
                         2.0 * p_i * sinh(beta*tau) * std::pow(beta,n) ;
          // If at boundary, cut in half b/c rectangles only half normal size
          return (beta == 0 or beta == (p.size()-1)*delta_b) ? 0.5*val : val;});

  return delta_b * ranges::accumulate(funcVal,0.0);
  // return the sum at all requested points, 
  // multiplied by the width of the rectangles
  // to give the Riemann summed area

}

#endif
