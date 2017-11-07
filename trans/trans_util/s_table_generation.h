#include <iostream>
#include <cmath>
#include <vector>
#include "bessel_K1.h" 


auto free_gas_s_table(const double& trans_weight, const double& alpha_sc, 
  const int& ndmax, const double& delta, std::vector<double>& s_free, 
  std::vector<double>& ap ){
  /* This creates a table of translational S(a,-b) values in the array s_free, 
   * where the table is evaluated at beta spaces of delta. This translational 
   * S(a,-b) can be used for free gasses, but can also represent the 
   * translational component for liquid moderators like water.  
   */
  double beta = 0, wal = trans_weight * alpha_sc;
  int j = 0;

  while ( true ){
    // Eq. 533 in the NJOY manual
    s_free[j] = exp( - 0.25 * (wal-beta) * (wal-beta) / wal ) / 
                  sqrt( 4 * M_PI * wal );
    ap[j] = beta;
    beta += delta; j += 1;
    // nsd must always be odd for use with Simpson's rule.
    if ( (j%2 == 1) and ( j+1 >= ndmax or 1e-7*s_free[0] >= s_free[j-1] ) ){
      return j;
    } // if break
  } // while
}


auto diffusion_s_table( const double& trans_weight, const double& alpha_sc, 
  const int& ndmax, const double& delta, std::vector<double>& s_diffusion, 
  std::vector<double>& ap, const double& diffusion){
  /* This creates a table of translational S(a,-b) values in the array s_free, 
   * where the table is evaluated at beta spaces of delta. This translational 
   * S(a,-b) can be used for free gasses, but can also represent the 
   * translational component for liquid moderators like water.  
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

    ap[j] = beta;
    beta += delta; j += 1;

    // nsd must always be odd for use with Simpson's rule.
    if (j%2 == 1 and ( (j+1 >= ndmax ) or 
                       (1e-7*s_diffusion[0]>=s_diffusion[j-1]) ) ){ 
      return j;
    } // if break
  } // while
}

