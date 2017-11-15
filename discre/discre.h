#include <iostream>
#include <vector>
#include <cmath>
#include "discre_util/bfill.h"
#include "discre_util/exts.h"
#include "discre_util/prepareParams.h"
#include "discre_util/oscLoopFuncs.h"

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
    
    int nn = oscillatorLoop( alpha, dbw, ar, scaling, wts, wtn, bes, ben, 
      energyNorm, a, maxdd, energy.size(), wt, tbart, weights, dist, 
      temp_vec[itemp] );

    // Sort the discrete lines, and throw out the smallest ones
    int n = nn;
    for ( auto i = 1; i < n; ++i ){
      for ( auto j = i+1; j < n+1; ++j ){
        if ( wts[j] > wts[i] ){
          auto save = wts[j];
          wts[j] = wts[i];
          wts[i] = save;
          save = bes[j];
          bes[j] = bes[i];
          bes[i] = save;
        } 
      }
    }
    
    int i = 0;
    int idone = 0;
    while ( i < nn ){
      i += 1;
      n = i;
      if ( wts[i-1] < 1.0e-6 and i > 5 ){ break; }
    }

    for ( auto m = 0; m < n; ++m ){
      for ( auto j = 0; j < beta.size(); ++j ){
        auto be = -betan[j] - bes[m];
       // std::cout << be << std::endl;
      } 
    }
    
    //for ( auto entry : bes ){ std::cout << entry << std::endl; }
    //std::cout << nn << std::endl;
    return;
  }


}
