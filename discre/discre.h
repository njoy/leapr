#include <iostream>
#include <vector>
#include "discre_util/bfill.h"

auto discre(const double& sc, const double& scaling, 
  const std::vector<double>& alpha, 
  const std::vector<double>& beta, const double& tev, const double& lambda_s, 
  std::vector<double>& osc_energies, std::vector<double>& osc_weights,
  const double& tbeta, std::vector<double>& t_eff_vec, 
  const std::vector<double>& temp_vec, int itemp ){
  int maxbb = 2 * beta.size() + 1;
  int maxdd = 500;
  double bk = 8.617385E-5;

  // Set up oscillator parameters
  int weight = 0;
  std::vector<double> osc_energies_norm ( osc_energies.size(), 0.0 );
  for ( auto i = 0; i < osc_energies.size(); ++i ){
    osc_energies_norm[i] = osc_energies[i] / tev;
    weight += osc_weights[i];
  }
  double tsave = 0.0;
  
  std::vector<double> eb(50, 0.0), ar(50, 0.0), dist(50, 0.0), dbw(50, 0.0);
  double sn, cn;
  for ( auto i = 0; i < osc_energies.size(); ++i ){
    eb[i] = exp( osc_energies_norm[i] / 2 );
    sn = (eb[i] - (1/eb[i]))/2;
    cn = (eb[i] + (1/eb[i]))/2;
    ar[i] = osc_weights[i] / ( sn * osc_energies_norm[i] );
    dist[i] = osc_weights[i] * osc_energies[i] * cn / ( 2*sn );
    tsave += dist[i] / bk;
    dbw[i] = ar[i] * cn;
  }
  
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
