#include "contin/contin_util/start_util/fsum.h"
#include <range/v3/all.hpp>

template <typename Float, typename Range>
auto normalize( Range beta_P, const Float& continWgt ){

  /* Rearranging Eq. 507 to get a definition for rho(beta), this is the 
   * equation that is being normalized to integrate to tbeta. 
   *
   *             rho(beta) = P(beta) * 2 * beta * sinh( beta / 2 )
   *
   * Inputs
   * ------------------------------------------------------------------------
   * p         : the vector to be normalized. This is P(beta), which is 
   *             defined by Eq. 507. 
   * delta_b   : spacing used in the Riemann sum for estimating the integral
   * continWgt : value to which the above equation should integrate
   *
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Integrate (using fsum) the above equation from 0 to infty.
   * * Scale each P(beta) value so that it is normlized to tbeta
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * P(beta) is amended
   */

  using std::sinh;

  auto integrand = [](auto xy){
    auto beta = std::get<0>(xy); 
    auto P = std::get<1>(xy); 
    return P*2.0*beta*sinh(beta*0.5);};

  Float invSum = continWgt / trapezoidIntegral(beta_P,integrand);
  return  ranges::view::values(beta_P)
        | ranges::view::transform([invSum](auto x){return x*invSum;});
}


