#include <iostream>
#include <vector>
#include <cmath>
#include "discre_util/bfill.h"
#include "discre_util/exts.h"
#include "discre_util/prepareParams.h"
#include "discre_util/oscLoopFuncs.h"
#include "discre_util/sint.h"
#include "discre_util/addDeltaFuncs.h"

void swap( double& a, double& b ){
  double c = a; a = b; b = c;
}

auto discre(const double& sc, const double& scaling, 
  const std::vector<double>& alpha, 
  const std::vector<double>& beta, const double& tev, const double& lambda_s, 
  std::vector<double>& energy, std::vector<double>& weights,
  const double& tbeta, std::vector<double>& t_eff_vec, 
  const std::vector<double>& temp_vec, int itemp,
  std::vector<std::vector<std::vector<double>>>& sym_sab,double twt  ){

  int maxbb = 2 * beta.size() + 1, maxdd = 500;
  double bk = 8.617385E-5;

  // Set up oscillator parameters
  // Prepare functions of beta
  double weight, tsave;

  std::vector<double> ar(50,0.0), t_eff_consts(50,0.0), lambda_i(50,0.0), 
    betaVals(50,0.0), exb(beta.size(),0.0), betan(beta.size(),0.0);

  prepareParams(energy, weights, tev, betaVals, weight, tsave, ar, t_eff_consts, lambda_i,
    bk, exb, betan, beta, sc );
  /* --> ar = [ weight / ( sinh( 0.5 * energy / tev ) * energy / tev ) ]
   *            This ends up being argument for bessel function in Eq. 537
   * --> betaVals = [ energy / tev ]
   * --> t_eff_consts = [ 0.5 * weight * energy / tanh( 0.5 * energy / tev ) ]
   *             This is used to calculate the effective temperature Eq. 544
   * --> lambda_i = [ weight / ( tanh( 0.5 * energy / tev ) * energy / tev ) ]
   *             This is lambda_i, defined in Eq. 538. Used for Eq. 537.
   */

  std::vector<double> bex( maxbb, 0.0 ), rdbex( maxbb, 0.0 );
  int nbx = bfill( bex, rdbex, betan );
  double wt = tbeta;
  double tbart = t_eff_vec[itemp]/temp_vec[itemp];
     
  // Main alpha loop
  for ( auto a = 0; a < alpha.size(); ++a ){

   
    // Get all sym_sab entries for a given alpha and temperature (vary beta)
    // for use in exts
    std::vector<double> input ( beta.size(), 0.0 );
    for ( auto b = 0; b < beta.size(); ++b ){
      input[b] = sym_sab[a][b][itemp];
    }

    std::vector<double> sex = exts( input, exb, betan );

    // Initialize delta loop
    std::vector<double> sexpb(beta.size(),0.0), bes(maxdd,0.0), wts(maxdd,0.0);
    
    int nn = oscillatorLoop( alpha, lambda_i, ar, scaling, wts, bes,  
      betaVals, a, maxdd, energy.size(), wt, tbart, weights, t_eff_consts, 
      temp_vec[itemp] );

    // Sort the discrete lines, and throw out the smallest ones
    int n = nn; double save;
    for ( auto i = 1; i < n; ++i ){
      for ( auto j = i+1; j < n+1; ++j ){
        if ( wts[j] > wts[i] ){
          swap( wts[j], wts[i] );
          swap( bes[j], bes[i] );
        } 
      }
    }
    
    int i = 0;
    while ( i < nn ){
      i += 1;
      n = i;
      if ( wts[i-1] < 1.0e-6 and i > 5 ){ break; }
    }

    // Add the continuous part to the scattering law
    for ( auto m = 0; m <= n; ++m ){
      for ( auto b = 0; b < beta.size(); ++b ){
        auto beta_val = -betan[b] - bes[m];
        auto add = wts[m] * sint(beta_val, bex, rdbex, sex, betan, b, 
            alpha[a], tbeta + twt, tbart, nbx);
        if ( add >= 1.0e-20 ){ sexpb[b] += add; }
      } 
    }

    // Add the delta functions to the scattering law
    double dwf = exp( -alpha[a]*scaling*lambda_s );
    addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n ); 

    // Record the results
    for ( auto j = 0; j < betan.size(); ++j ){
      sym_sab[a][j][itemp] = sexpb[j];
    }
  }
}
