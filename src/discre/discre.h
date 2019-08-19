#include "discre_util/bfill.h"
#include "discre_util/exts.h"
#include "discre_util/prepareParams.h"
#include "discre_util/oscLoopFuncs.h"
#include "discre_util/sint.h"
#include "discre_util/addDeltaFuncs.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>
#include <iostream>

template <typename Float>
void swap( Float& a, Float& b ){ Float c = a; a = b; b = c; }

template <typename Float, typename Range, typename RangeZipped>
auto discre( const Float& lambda_s, 
  const Float& twt, const Float& tbeta, Range alpha, Range beta, 
  const Float& temp, RangeZipped oscEnergiesWeights, Float& t_eff, Range& sab ){

  //for ( auto& a : alpha ){ a *= scaling; }
  //for ( auto& b : beta  ){ b *= sc;      }



  int maxdd = 500;

  // Set up oscillator parameters
  Range ar(50), t_eff_consts(50), debyeWaller(50), oscBetas(50), exb(beta.size());
  Float tev = temp*kb;

  prepareParams(oscEnergiesWeights, tev, oscBetas, ar, t_eff_consts,
    debyeWaller, exb, beta );

  /* --> ar = [ weight / ( sinh( 0.5 * energy / tev ) * energy / tev ) ]
   *            This ends up being argument for bessel function in Eq. 537
   * --> t_eff_consts = [ 0.5 * weight * energy / tanh( 0.5 * energy / tev ) ]
   *             This is used to calculate the effective temperature Eq. 544
   * --> debyeWaller = [ weight / ( tanh(0.5 * energy / tev) * energy / tev ) ]
   *             This is lambda_i, defined in Eq. 538. Used for Eq. 537.
   * --> exb = [ exp( -beta * sc / 2 ) ]
   *          This is used in calculating the sex vector, since to go from 
   *          S(a,b) --> S(a,-b) you need to multiply by exp( -beta )
   */

  Range rdbex( 2 * beta.size() + 1 );
  auto output = bfill(rdbex, beta);
  int nbx   = std::get<0>(output);
  Range bex = std::get<1>(output);

  Float tbart = t_eff/temp;
     
  for ( size_t a = 0; a < alpha.size(); ++a ){
    // Get all sab entries for a given alpha and temperature (vary beta)
    Range input ( beta.size() );
    for ( size_t b = 0; b < beta.size(); ++b ){
      input[b] = sab[b+a*beta.size()];
    }
    // input = sab values for constant temp and alpha. 
    // exb   = exp( -beta * sc / 2 ), which (following Eq. 509) we need in 
    //         order to go between S(a,b) and S(a,-b) 
    Range sex = exts( input, exb, beta );
    // sex is populated with sab entries, such that 
    //        sex = [ s3 s2 s1 s1 s2*exp(-beta) s3*exp(-beta) 0 ]
    //                             or 
    //         sex = [ s3 s2 s1 s2*exp(-beta) s3*exp(-beta) 0 0 ]
    //                (dependinng on first beta value)
    // The exp(-beta) values are explained above, because of Eq. 509

    // Initialize delta loop
    Range bes(maxdd,0.0), wts(maxdd,0.0);
    
    unsigned int nn = oscillatorLoop( alpha[a], debyeWaller, ar, wts, bes,  
      oscBetas, tbart, t_eff_consts, temp );

    // oscillator loop is mean to, for a given alpha and beta, populate the wts
    // vector with entries of W_k(alpha) for various k (see Eq. 542) and to
    // populate bes with entries of beta_k again for various k.
    // So at this point all we need to do is sum over k for W_k*S(a,b-b_k)


    // Sort the discrete lines, and throw out the smallest ones
    // Except for the first value, we're sorting wts and bes so that wts values
    // are in decreasing order.
    unsigned int n = nn; 

    for ( size_t i = 1; i < nn; ++i ){
      for ( size_t j = i+1; j <= n; ++j ){
        if ( wts[j] >= wts[i] ){
          swap( wts[j], wts[i] );
          swap( bes[j], bes[i] );
        } 
      }
    }
    for (size_t i = 1; i <= nn; ++i){
      n = i;
      if (wts[i-1] < 1e-6 and i > 5) break; 
    }

    // Add the continuous part to the scattering law
    Range sexpb(beta.size(),0.0);

    for ( size_t m = 0; m < n; ++m ){
      for ( size_t b = 0; b < beta.size(); ++b ){
        //auto beta_val = -beta[b] - bes[m];
        // This is explicitly evaluating Eq. 542, where wts is W_k(alpha), and
        // bes is a vector populated with beta_k values. sint is used to 
        // interpolate for the beta - beta_k piece of the equation.
        auto add = wts[m] * sint(-beta[b] - bes[m], bex, rdbex, sex, beta, b, 
            alpha[a]*(tbeta + twt), tbart, nbx);
        if ( add >= 1.0e-20 ){ sexpb[b] += add; }
      } 
    }

    // Add the delta functions to the scattering law
    Float debyeWallerExponential = exp( -alpha[a]*lambda_s );
    addDeltaFuncs( twt, debyeWallerExponential, bes, beta, wts, sexpb, n ); 

    // Record the results
    for ( size_t b = 0; b < beta.size(); ++b ){
      sab[b+a*beta.size()] = sexpb[b];
    }
  }
  double tsave = 0.0;
  for ( size_t n = 0; n < t_eff_consts.size(); ++n){
    tsave += t_eff_consts[n]/kb;
  }

  t_eff = (tbeta+twt)*t_eff + tsave;
}














