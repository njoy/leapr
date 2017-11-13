#include <iostream>
#include <vector>
#include <cmath>
#include "discre_util/bfill.h"
#include "discre_util/exts.h"
#include "discre_util/prepareParams.h"
#include "discre_util/bfact.h"

auto discre(const double& sc, const double& scaling, 
  const std::vector<double>& alpha, 
  const std::vector<double>& beta, const double& tev, const double& lambda_s, 
  std::vector<double>& energy, std::vector<double>& weights,
  const double& tbeta, std::vector<double>& t_eff_vec, 
  const std::vector<double>& temp_vec, int itemp,
  std::vector<std::vector<std::vector<double>>>& sym_sab ){

  int maxbb = 2 * beta.size() + 1, maxdd = 500;
  double bk = 8.617385E-5;

  // Set up oscillator parameters
  // Prepare functions of beta
  double weight, tsave;
  std::vector<double> ar(50,0.0), dist(50,0.0), dbw(50,0.0), 
    energyNorm(energy.size(),0.0), exb(beta.size(),0.0), betan(beta.size(),0.0);

  prepareParams(energy, weights, tev, energyNorm, weight, tsave, ar, dist,dbw,
    bk, exb, betan, beta, sc );

  std::vector<double> bex( maxbb, 0.0 ), rdbex( maxbb, 0.0 );
  auto ndx = bfill( bex, rdbex, betan );
  double wt = tbeta;
  double tbart = t_eff_vec[itemp]/temp_vec[itemp];
     
  // Main alpha loop
  for ( auto a = 0; a < alpha.size(); ++a ){
    double dwf = exp( -alpha[a]*scaling*lambda_s );
    std::vector<double> sex ( maxbb, 0.0 );
    std::vector<double> input ( beta.size(), 0.0 );
    for ( auto b = 0; b < beta.size(); ++b ){
      input[b] = sym_sab[a][b][itemp];
    }
    exts( input, sex, exb, betan );

    // Initialize delta loop
    std::vector<double> sexpb (beta.size(), 0.0);
    std::vector<double> ben ( maxdd, 0.0 );
    std::vector<double> bes ( maxdd, 0.0 );
    std::vector<double> wtn ( maxdd, 0.0 );
    std::vector<double> wts ( maxdd, 0.0 );
    wtn[0] = 1;
    int nn = 1;
    int n = 0;
    
    // Loop over all oscillators
    for ( auto i = 0; i < energy.size(); ++i ){
      double dwc = alpha[a]*scaling*dbw[i];
      double x   = alpha[a]*scaling*ar[i];
      std::vector<double> bminus (50,0.0), bplus(50,0.0);
      double bzero = bfact( x, dwc, energy[i], bplus, bminus );

      // do convolution for delta function
      for ( auto m = 0; m < nn; ++m ){
        double besn = ben[m];
        double wtsn = wtn[m]*bzero;
        if ( besn <= 0 or wtsn >= 1e-8 ){
          if ( n < maxdd ){
            bes[n] = besn;
            wts[n] = wtsn;
            n += 1;
          }
        }
      }
      
      // negative n terms


      return;

    }   
    return;
  }


}
