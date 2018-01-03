#include <iostream> 
#include <cmath>
#include <vector>
#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"


void trans( const std::vector<double>& alpha, const std::vector<double>& beta,
  const double& trans_weight, double delta, const double& diffusion, 
  const double& sc, const double& scaling, const int& itemp, 
  const double& lambda_s, const double& tbeta, std::vector<double>& t_eff_vec, 
  const std::vector<double>& temp_vec, 
  std::vector<std::vector<std::vector<double>>>& sym_sab ){

  /* Overview
   * ------------------------------------------------------------------------
   * trans is used to control the addition of a translational contribution to
   * the presently continuous S(a,b). This translational component can be 
   * either diffusive for a liquid or a free gas. This contribution is 
   * calculated using s_table_generation, the existant beta grid is shifted 
   * to match the contribution, and then the two pieces are convolved with the
   * result recorded. 
   *
   * Inputs
   * ------------------------------------------------------------------------
   * alpha        : vector of alpha values
   * beta         : vector of beta values
   * trans_weight : translational weight, appears as w_t in the translational
   *                equations (pg. 654), and is user provided from Card13.
   * delta        : user-defined spacing. Originally provided as a spacing in
   *                eV, but was made nondimensional in leapr.cpp. Corresponds
   *                to preliminary beta grid spacing.
   * diffusion    : diffusion contant (zero for free gas). Appears as c in the
   *                translational equations, and is user provided from Card13.
   * sc           : scaling parameter used for beta, defined in leapr.cpp
   * scaling      : scaling parameter used for alpha, defined in leapr.cpp
   * itemp        : current temp index
   * lambda_s     : current debye-waller coefficient
   * tbeta        : normalization for the continuous part
   * t_eff_vec    : effective temperature vector
   * temp_vec     : temperature vector
   * sym_sab      : current scattering law, which at this point only has the
   *                solid / continuous freq. distribution contribution
   * 
   * Operations
   * ------------------------------------------------------------------------
   *
   *  
   * Outputs
   * ------------------------------------------------------------------------
   */



  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;

  std::vector<double> sabTrans(ndmax,0.0), ap(ndmax,0.0), sab(ndmax, 0.0),
    betan(beta.size(),0.0);

  double nsd, alpha_sc, ded;
  // loop over alpha values
  for ( auto a = 0; a < alpha.size(); ++a ){
    alpha_sc = alpha[a] * scaling;

    ded = diffusion == 0 ? 
      0.2 * sqrt( trans_weight * alpha_sc ) :
      0.4 * trans_weight * diffusion * alpha_sc / 
        sqrt( 1.0 + 1.42*trans_weight*diffusion*diffusion*alpha_sc );

    delta = std::min( ded, 10.0 * alpha_sc * delta );

    nsd = diffusion == 0 ? free_gas_s_table ( trans_weight, alpha_sc, ndmax, 
                                             delta, sabTrans ) : 
                           diffusion_s_table( trans_weight, alpha_sc, ndmax, 
                                             delta, sabTrans, diffusion );
    if ( nsd > 1 ){
      for ( int b = 0; b < beta.size(); ++b ){
        betan[b] = beta[b] * sc;
        ap[b] = sym_sab[a][b][itemp];
      }

      // loop over beta values
      for ( int b = 0; b < beta.size(); ++b ){
        double be = betan[b];

        // prepare table of continuous ss on new interval
        sbfill( sab, nsd, delta, be, ap, betan, ndmax );

        // convolve s-transport with s-continuous
        double s = 0;
        for ( int i = 0; i < nsd; ++i ){
          // f = 2 if i is even, 4 if i is odd, except 1 at boundaries
          double f = 2*(i%2)+2;
          if ( i == 0 or i == nsd - 1 ){ f = 1; }
                    
          s = s + f*sabTrans[i]*sab[nsd+i-1] + f*sabTrans[i]*sab[nsd-i-1]*exp(-i*delta);

        }
        s = s < 1e-30 ? 0 : s * delta / 3;
                
        double st = terps(sabTrans,nsd,delta,be);

	// This accounts for the first term in Eq. 535
        if ( st > 0.0 ){ s += exp( -alpha_sc * lambda_s ) * st; }

        sym_sab[a][b][itemp] = s;

      } // for beta
    } // if nsd > 0
  } // for alpha
  
  // Update the effective temperature, following Eq. 536
  t_eff_vec[itemp] = (tbeta*t_eff_vec[itemp] + trans_weight*temp_vec[itemp]) /
                     ( tbeta + trans_weight );


}



