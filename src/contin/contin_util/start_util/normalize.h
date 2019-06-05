#include "contin/contin_util/start_util/fsum.h"
#include <range/v3/all.hpp>

template <typename Float, typename Range>
auto normalize( Range P_beta, const Float& continWgt ){

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

  Float invSum = continWgt / fsum( 1, P_beta, 0.5 );  
  return  ranges::view::values(P_beta)
        | ranges::view::transform([invSum](auto x){return x*invSum;});
}


