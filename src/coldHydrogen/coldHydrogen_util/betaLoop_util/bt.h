#include <range/v3/all.hpp>

template <typename Float>
auto bt ( int j, const Float& x ){
  /* Inputs
   * ------------------------------------------------------------------------
   * j  : This calculates the jth statistical weight factor
   * x  : Originally from coldh.h, it's set equal to de / tev, which is either
   *      0.0147 / kb*T or 0.0074 / kb*T depending on whether its hydrogen
   *      or deuterium, respectively.
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * This calculates the statistical weight factor Pj for cold hydrogen or 
   * deuterium calculation. This goes in Eq. 567 - 568, and is called by the
   * betaLoop function.
   * 
   *                    ( 2j + 1 ) * exp( - j*(j+1)*x/2 ) 
   *   Pj = -------------------------------------------------------------------
   *        2( (2*1+1)e^(-1*2*x) + (2*3+1)e^(-3*4*x) + (2*5+1)e^(-5*6*x) + ... )
   *
   * if j is odd ( ortho hydrogen or para deuterium )
   *
   *                    ( 2j + 1 ) * exp( - j*(j+1)*x/2 ) 
   *   Pj = -------------------------------------------------------------------
   *        2( (2*2+1)e^(-2*3*x) + (2*4+1)e^(-4*5*x) + (2*6+1)e^(-6*7*x) + ... )
   * if j is even ( para hydrogen or ortho deuterium )
   *
   */

  // Loop either 0 --> 18 or 1 --> 19 to calculate the denominator first
  auto kRange = ranges::view::iota(0,10) 
              | ranges::view::transform([j](auto i){ 
                  return j%2 == 0 ? i*2 : i*2 + 1; });
  auto b = ranges::accumulate( kRange | ranges::view::transform( [x](auto k){ 
                           return (2*k+1) * exp( -k * (k+1) * x / 2 ); } ),0.0);

  return (j+0.5) * exp(-j*(j+1)*x*0.5) / b;
}

