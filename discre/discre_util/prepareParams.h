#include <iostream>
#include <vector>
#include <cmath>


auto prepareParams( const std::vector<double>& energy, 
  const std::vector<double>& weights, const double& tev, 
  std::vector<double>& betaVals, double& weight, 
  double& tsave, std::vector<double>& ar, std::vector<double>& dist, 
  std::vector<double>& dbw, const double& bk, 
  std::vector<double>& exb, std::vector<double>& betan, 
  const std::vector<double>& beta, const double& sc ){

  // Set up oscillator parameters
  
  weight = 0.0;
  tsave = 0.0;
  for ( auto i = 0; i < energy.size(); ++i ){
    betaVals[i] = energy[i] / tev;
    weight += weights[i];

    ar[i]   = weights[i] / ( sinh(0.5*betaVals[i]) * betaVals[i] );
    dist[i] = 0.5 * weights[i] * energy[i] / tanh(0.5*betaVals[i]);
    dbw[i]  = weights[i] / ( tanh(0.5*betaVals[i]) * betaVals[i] );
    tsave  += dist[i] / bk;
  }

  for ( auto b = 0; b < betan.size(); ++b ){
    exb[b] = exp( -beta[b]*sc*0.5 );
    betan[b] = beta[b]*sc;
  } 
}
