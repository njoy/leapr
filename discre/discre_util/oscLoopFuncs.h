#include <iostream>
#include <vector>
#include "bfact.h"

void posNegTerms( int i, int& n, const double& normalizedEnergy, 
  const std::vector<double>& b_minus_or_plus,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, const int nn,
  int pos_or_neg ){

  int k = 0;

  while ( k < 50 ){
    k += 1;
    if ( b_minus_or_plus[k-1] <= 0 ){ return; } 
    for ( auto m = 0; m <= nn; ++m ){
      if ( wtn[m] * b_minus_or_plus[k-1] >= 1e-8 and n < bes.size() ){
        n += 1;
        bes[n] = ben[m] + pos_or_neg * k * normalizedEnergy;
        wts[n] = wtn[m] * b_minus_or_plus[k-1];
      }
    }
  }
}


auto oscillatorLoop( const std::vector<double>& alpha,
  std::vector<double>& dbw, std::vector<double>& ar, const double& scaling,
  std::vector<double>& wts, std::vector<double>& wtn, std::vector<double>& bes,
  std::vector<double>& ben, std::vector<double>& energyNorm, int a, int maxdd,
  int numOscillators, double& wt, double& tbart, std::vector<double>& weight,
  std::vector<double>& dist, const double& temp ){
  double bk = 8.617385E-5;
  int nn = 1;
  int n;
  // Loop over all oscillators
  for ( auto i = 0; i < numOscillators; ++i ){
    n = 0;
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
    posNegTerms( i, n, energyNorm[i], bminus, wts, wtn, bes, ben, nn, -1 );

    // positive n terms
    posNegTerms( i, n, energyNorm[i], bplus, wts, wtn, bes, ben, nn, 1 );

    // continue loop
    nn = n;
    for ( auto m = 0; m <= nn; ++m ){
      ben[m] = bes[m];
      wtn[m] = wts[m];
    }
    wt += weight[i];
    tbart += dist[i] / ( bk * temp );

  }   
  return nn;
}


