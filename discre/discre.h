#include <iostream>
#include <vector>
#include <cmath>
#include "discre_util/bfill.h"
#include "discre_util/exts.h"
#include "discre_util/prepareParams.h"

auto discre(const double& sc, const double& scaling, 
  const std::vector<double>& alpha, 
  const std::vector<double>& beta, const double& tev, const double& lambda_s, 
  std::vector<double>& energy, std::vector<double>& weights,
  const double& tbeta, std::vector<double>& t_eff_vec, 
  const std::vector<double>& temp_vec, int itemp ){
  int maxbb = 2 * beta.size() + 1;
  int maxdd = 500;
  double bk = 8.617385E-5;

  // Set up oscillator parameters
  double weight, tsave;
  std::vector<double> energyNorm ( energy.size(), 0.0 );
  std::vector<double> ar(50, 0.0), dist(50, 0.0), dbw(50, 0.0);
  prepareParams(energy, weights, tev, energyNorm, weight, tsave,ar, dist,dbw, bk );
//  for ( auto i = 0; i < energy.size(); ++i ){
//    energyNorm[i] = energy[i] / tev;
//    weight += weights[i];
//  }

//  double tsave = 0.0;
  
//  double sinh_term, cosh_term, tanh_term;
//  for ( auto i = 0; i < energy.size(); ++i ){
//    sinh_term = sinh( 0.5*energyNorm[i] );
//    cosh_term = cosh( 0.5*energyNorm[i] );
//    tanh_term = tanh( 0.5*energyNorm[i] );
//    ar[i]   = weights[i] / ( sinh_term * energyNorm[i] );
//    dist[i] = 0.5 * weights[i] * energy[i] / tanh_term;
//    tsave += dist[i] / bk;
//    dbw[i] = ar[i] * cosh_term;
//  }
  
  // Prepare functions of beta
  std::vector<double> exb ( beta.size(), 0.0 ), betan ( beta.size(), 0.0 );
  for ( auto i = 0; i < beta.size(); ++i ){
    exb[i] = exp(-beta[i]*sc/2); 
    betan[i] = beta[i]*sc;
  }
  std::vector<double> bex( maxbb, 0.0 ), rdbex( maxbb, 0.0 );
  auto ndx = bfill( bex, rdbex, betan );
  double wt = tbeta;
  double tbart = t_eff_vec[itemp]/temp_vec[itemp];
     
  // Main alpha loop
  for ( auto a = 0; a < alpha.size(); ++a ){
    double dwf = exp( -alpha[a]*scaling*lambda_s );
    std::cout << dwf << std::endl;
  }

}
