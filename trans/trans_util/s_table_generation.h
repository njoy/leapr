#include <iostream>
#include <cmath>
#include <vector>
#include "bessel_K1.h" 


auto free_gas_s_table(const double& trans_weight, const double& alpha_sc, 
  const int& ndmax, const double& delta, std::vector<double>& s_free ){
  /* Overview
   * ------------------------------------------------------------------------
   * This creates a table of translational S(a,-b) values in the array s_free, 
   * where the table is evaluated at beta spaces of delta. This translational 
   * S(a,-b) can be used for free gasses, but can also represent the 
   * translational component for liquid moderators like water.  
   *
   * Inputs
   * ------------------------------------------------------------------------
   * trans_weight : translational weight, which appears as w_t in equations 
   *                in the Diffusion and Free-Gas Translation part of the 
   *                manual (pg. 654). User provided value from Card13.
   * alpha_sc     : scaled alpha value
   * ndmax        : maximum number of values in the sd matrix 
   * delta        : ideal delta value computed for spacing S_t(a,-b)
   * s_free       : empty vector (S_t(a,-b))
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Calculate the right hand side of Eq. 533 for increasing values of beta,
   *   and insert them in s_free
   * * Exit loop when either s_free is full or when the first and most recent
   *   values of beta are different enough
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * the s_free vector is filled with values corresponding to Eq. 533. 
   *   Notice that these values are for S_t(a,-b)
   */


  double beta = 0, wal = trans_weight * alpha_sc;
  int j = 0;

  while ( true ){
    // Eq. 533 in the NJOY manual
    s_free[j] = exp( - std::pow(wal - beta, 2) / (4 * wal ) ) / 
                 sqrt( 4 * M_PI * wal );
    beta += delta; j += 1;
    // nsd must always be odd for use with Simpson's rule.
    if ( (j%2 == 1) and ( j+1 >= ndmax or 1e-7*s_free[0] >= s_free[j-1] ) ){
      return j;
    } // if break
  } // while
}


auto diffusion_s_table( const double& trans_weight, const double& alpha_sc, 
  const int& ndmax, const double& delta, std::vector<double>& s_diffusion, 
  const double& diffusion){
  /* Overview
   * ------------------------------------------------------------------------
   * This creates a table of translational S(a,-b) values in the array s_free, 
   * where the table is evaluated at beta spaces of delta. This translational 
   * S(a,-b) can be used for free gasses, but can also represent the 
   * translational component for liquid moderators like water.  
   *
   * Inputs
   * ------------------------------------------------------------------------
   * trans_weight : translational weight, which appears as w_t in equations 
   *                in the Diffusion and Free-Gas Translation part of the 
   *                manual (pg. 654). User provided value from Card13.
   * alpha_sc     : scaled alpha value
   * ndmax        : maximum number of values in the sd matrix 
   * delta        : ideal delta value computed for spacing S_t(a,b)
   * s_diffusion  : empty vector (S_t(a,b))
   * diffusion    : diffusion constant, which appears in equations as c, and
   *                if a user provided value from Card13.
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Calculate the right hand side of Eq. 531 for increasing values of beta,
   *   and insert them in s_diffusion
   * * Exit loop when either s_diffusion is full or when the first and most 
   *   recent values of beta are different enough
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * the s_diffusion vector is filled with values corresponding to Eq. 531.
   */

  double beta = 0;
  int j = 0;
  double wda = trans_weight * diffusion * alpha_sc; 
  while ( true ){
    
    double expTerm = sqrt(beta*beta + 4*wda*wda) * 
                     sqrt( diffusion*diffusion + 0.25 );

    // This is evaluating Eq. 531 from the NJOY manual        
    s_diffusion[j] = ( 2*wda/M_PI ) * sqrt( diffusion*diffusion + 0.25 ) * 
                     bessel_K1_gen( expTerm ) / sqrt( beta*beta + 4*wda*wda );
    s_diffusion[j] *= expTerm <= 1 ? exp( 2*wda*diffusion + beta/2 ) :
                                     exp( 2*wda*diffusion + beta/2 - expTerm); 

    beta += delta; j += 1;

    // nsd must always be odd for use with Simpson's rule.
    if ( ( (j+1 >= ndmax) or (1e-7*s_diffusion[0] >= s_diffusion[j-1]) ) and
       ( j%2 == 1 ) ){
      return j;
    } // if break
  } // while
}

