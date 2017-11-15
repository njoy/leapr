#include <iostream>
#include <vector>
#include "bfact.h"

void negativeTerms( int i, int& n, const double& normalizedEnergy, 
  const std::vector<double>& bminus,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, const int nn){

  int k = 0;

  while ( k < 50 ){
    k += 1;
    if ( bminus[k-1] <= 0 ){ return; } 
    for ( auto m = 0; m <= nn; ++m ){
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
    for ( auto m = 0; m <= nn; ++m ){
      if ( wtn[m] * bplus[k-1] >= 1e-8 and n < bes.size() ){
        n += 1;
        bes[n] = ben[m] + k * normalizedEnergy;
        wts[n] = wtn[m] * bplus[k-1];
      }
    }
  }
}

void oscillatorLoop( std::vector<double>& alpha,
  std::vector<double>& dbw, std::vector<double>& ar, const double& scaling,
  std::vector<double>& wts, std::vector<double>& wtn, std::vector<double>& bes,
  std::vector<double>& ben, std::vector<double>& energyNorm, int a, int maxdd,
  int numOscillators, double& wt, double& tbart, std::vector<double>& weight,
  std::vector<double>& dist, double& temp ){
  double bk = 8.617385E-5;
  int nn = 1;

  // Loop over all oscillators
  for ( auto i = 0; i < numOscillators; ++i ){
    int n = 0;
    double dwc = alpha[a]*scaling*dbw[i];
    double x   = alpha[a]*scaling*ar[i];
    std::vector<double> bminus (50,0.0), bplus(50,0.0);
    double bzero = bfact( x, dwc, energyNorm[i], bplus, bminus );
    
    // do convolution for delta function
    for ( auto m = 0; m <= nn; ++m ){
      double besn = ben[m];
      double wtsn = wtn[m]*bzero;
      if ( besn < 0 or wtsn >= 1e-8 ){
        if ( n < maxdd ){
         bes[n] = besn;
          wts[n] = wtsn;
          n += 1;
        }
      }
    }
    n -= 1;

    // negative n terms
    negativeTerms( i, n, energyNorm[i], bminus, wts, wtn, bes, ben, nn );

    // positive n terms
    positiveTerms( i, n, energyNorm[i], bplus, wts, wtn, bes, ben, nn );

    // continue loop
    nn = n;
    for ( auto m = 0; m <= nn; ++m ){
      ben[m] = bes[m];
      wtn[m] = wts[m];
    }
    wt += weight[i];
    tbart += dist[i] / ( bk * temp );


  }   
  return;
}


