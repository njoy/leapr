#include "contin/contin_util/start_util/fsum.h"

template <typename F, typename A>
void normalize( A& p, const F& tbeta, const A& betaGrid ){

  /* Rearranging Eq. 507 to get a definition for rho(beta), this is the 
   * equation that is being normalized to integrate to tbeta. 
   *
   *             rho(beta) = P(beta) * 2 * beta * sinh( beta / 2 )
   *
   * Inputs
   * ------------------------------------------------------------------------
   * p       : the vector to be normalized. This is P(beta), which is 
   *           defined by Eq. 507. 
   * delta_b : spacing used in the Riemann sum for estimating the integral
   * tbeta   : value to which the above equation should integrate
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

  F inv_sum = tbeta / fsum( 1, p, 0.5, betaGrid );  // This gives you 
  for ( F& entry : p ){ entry = entry * inv_sum; } 
}


