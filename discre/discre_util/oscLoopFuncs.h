#include <iostream>
#include <vector>


void negativeTerms( int i, int& n, const double& normalizedEnergy, 
  const std::vector<double>& bminus,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, const int nn){

  int k = 0;

  while ( k < 50 ){
    k += 1;
    if ( bminus[k-1] <= 0 ){ return; } 
    for ( auto m = 0; m < nn; ++m ){
      if ( wtn[m] * bminus[k-1] >= 1e-8 and n < bes.size() ){
        n += 1;
        bes[n] = ben[m] - k * normalizedEnergy;
        wts[n] = wtn[m] * bminus[k-1];
      }
    }
  }
}

void positiveTerms( int i, int& n, const double& normalizedEnergy, 
  const std::vector<double>& bplus,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, const int nn){

  int k = 0;

  while ( k < 50 ){
    k += 1;
    if ( bplus[k-1] <= 0 ){ return; } 
    for ( auto m = 0; m < nn; ++m ){
      if ( wtn[m] * bplus[k-1] >= 1e-8 and n < bes.size() ){
        n += 1;
        bes[n] = ben[m] + k * normalizedEnergy;
        wts[n] = wtn[m] * bplus[k-1];
      }
    }
  }
}
