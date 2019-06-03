#include "simple/generalTools/values.h"

template <typename V, typename F>
auto prepareParams( 
  const std::vector<std::tuple<double,double>>& oscEnergiesWeights,
  const F& tev, V& betaVals, F& weight, F& tsave, V& ar, V& dist, V& dbw, 
  V& exb, V& betan, const V& beta, const F& sc ){
  /* --> ar = [ weight / ( sinh( 0.5 * energy / kbT ) * energy / kbT ) ]
   *            This ends up being argument for bessel function in Eq. 537
   * --> betaVals = [ energy / kbT ]
   * --> t_eff_consts = [ 0.5 * weight * energy / tanh( 0.5 * energy / kbT ) ]
   *             This is used to calculate the effective temperature Eq. 544
   * --> lambda_i = [ weight / ( tanh( 0.5 * energy / kbT ) * energy / kbT ) ]
   *             This is lambda_i, defined in Eq. 538. Used for Eq. 537.
   * --> exb = [ exp( -beta * sc / 2 ) ]
   *          This is used in calculating the sex vector, since to go from 
   *          S(a,b) --> S(a,-b) you need to multiply by exp( -beta )
   */


  // Set up oscillator parameters
  
  weight = 0.0;
  tsave = 0.0;
  for ( size_t i = 0; i < oscEnergiesWeights.size(); ++i ){
    double E = std::get<0>(oscEnergiesWeights[i]);
    double w = std::get<1>(oscEnergiesWeights[i]);

    betaVals[i] = E / tev;
    weight += w;

    ar[i]   = w / ( sinh(0.5*betaVals[i]) * betaVals[i] );
    dist[i] = 0.5 * w * E / tanh(0.5*betaVals[i]);
    dbw[i]  = w / ( tanh(0.5*betaVals[i]) * betaVals[i] );
    tsave  += dist[i] / kb;
  }

  for ( size_t b = 0; b < betan.size(); ++b ){
    exb[b] = exp( -beta[b]*sc );
    betan[b] = beta[b]*sc;
  } 
}
