#include <iostream>
#include <vector>
#include <cmath>
#include "oscLoopFuncs_util/bfact.h"

void posNegTerms( int& n, const double& beta_i, 
  const std::vector<double>& b_minus_or_plus,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, const int nn,
  int pos_or_neg ){

  int k = 0;
  /* There are 50 entries in bplus and bminus (not all of them necessarily
   * nonzero). Loop through them. 
   */
  for ( auto k = 0; k < 50; ++k ){
    if ( b_minus_or_plus[k] <= 0 ){ return; } 
    for ( auto m = 0; m < nn; ++m ){
      if ( wtn[m] * b_minus_or_plus[k] >= 1e-8 and n < bes.size() ){
        n += 1;
        bes[n] = ben[m] + pos_or_neg * (k+1) * beta_i;
        wts[n] = wtn[m] * b_minus_or_plus[k];
      }
    }
  }
}


auto oscillatorLoop( const std::vector<double>& alpha,
  std::vector<double>& lambda_i, std::vector<double>& ar, const double& scaling,
  std::vector<double>& wts, std::vector<double>& bes,
  std::vector<double>& betaVals, int a, int maxdd,
  int numOscillators, double& wt, double& tbart, std::vector<double>& weight,
  std::vector<double>& dist, const double& temp ){
  /* Note that the betaVals vector is a vector whose entries are the oscillator
   * energies scaled by temp in eV, thus making them beta values.
   */

  std::vector<double> ben( maxdd, 0.0 ), wtn( maxdd, 0.0 ); wtn[0] = 1.0;

  double bk = 8.617385E-5, alpha_lambda_i, x, bzero;
  int n = 0, nn = 0;

  // Loop over all oscillators
  for ( auto i = 0; i < numOscillators; ++i ){
    nn = n + 1;

    alpha_lambda_i = alpha[a]*scaling*lambda_i[i];
    //             = alpha*scaling*weight / (tanh(0.5*energy/tev) * energy/tev)
    //             = scaled alpha * lambda_i (lambda_i defined in Eq. 538)
    
    x              = alpha[a]*scaling*ar[i];
    //             = alpha*scaling*weight / (sinh(0.5*energy/tev) * energy/tev)
    //             = argument for bessel function in Eq. 537

    /* bfact populates bplus and bminus with A_in terms from Eq. 537.
     * The nth entry of bplus or bminus corresponds to a specific alpha and 
     * i value. The bzero output is either 
     *                 I0(x)*e^(-alpha*lambda_i)
     * or 
     *                I0(x)*e^(-alpha*lambda_i+x)
     * depending on the size of x.
     */
    std::vector<double> bminus (50,0.0), bplus(50,0.0);
    bzero = bfact( x, alpha_lambda_i, betaVals[i], bplus, bminus );
    
    // do convolution for delta function
    n = 0;
    for ( auto m = 0; m < nn; ++m ){
      if ( (ben[m] <= 0 or wtn[m]*bzero >= 1e-8) and n < maxdd ){
        bes[n] = ben[m];
        wts[n] = wtn[m]*bzero;
        n += 1;
      }
    }
    n -= 1;

    // negative n terms
    posNegTerms( n, betaVals[i], bminus, wts, wtn, bes, ben, nn, -1 );

    // positive n terms
    posNegTerms( n, betaVals[i], bplus,  wts, wtn, bes, ben, nn, 1  );

    // continue loop
    // Copy first n entries of permanent array into our temporary arrays
    for ( auto m = 0; m <= n; ++m ){
      ben[m] = bes[m];
      wtn[m] = wts[m];
    }
    wt += weight[i];
    // Effective temperature is amended, this ( kind of ) follows Eq. 544.
    tbart += dist[i] / ( bk * temp );

  }   
  return n;
}


