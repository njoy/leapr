#include <iostream>
#include <vector>
#include <cmath>


auto prepareParams( std::vector<double>& energy, std::vector<double>& weights,
  const double& tev, std::vector<double>& energyNorm, double& weight, 
  double& tsave, std::vector<double>& ar, std::vector<double>& dist, 
  std::vector<double>& dbw, const double& bk ){

  // Set up oscillator parameters
  weight = 0.0;
  for ( auto i = 0; i < energy.size(); ++i ){
    energyNorm[i] = energy[i] / tev;
    weight += weights[i];
  }

  tsave = 0.0;
  
  double sinh_term, cosh_term, tanh_term;
  for ( auto i = 0; i < energy.size(); ++i ){
    sinh_term = sinh( 0.5*energyNorm[i] );
    cosh_term = cosh( 0.5*energyNorm[i] );
    tanh_term = tanh( 0.5*energyNorm[i] );
    ar[i]   = weights[i] / ( sinh_term * energyNorm[i] );
    dist[i] = 0.5 * weights[i] * energy[i] / tanh_term;
    tsave += dist[i] / bk;
    dbw[i] = ar[i] * cosh_term;
  }
 
  return ;
}
