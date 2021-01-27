#include "oscLoopFuncs_util/bfact.h"
#include <range/v3/all.hpp>
#include <iostream>
#include "generalTools/constants.h"

void posNegTerms( int& n, const double& beta_i, 
  const std::vector<double>& b_minus_or_plus, std::vector<double>& wts, 
  const std::vector<double>& wtn, std::vector<double>& bes, 
  const std::vector<double>& ben, const int nn, int pos_or_neg ){

  /* The point of these loops are to prepare wts and bes for the summation in
   * Eq. 542, which involves a sum of sums of sums etc, so we're really 
   * interested in getting cross terms here. 
   *
   * For the first oscillator (i=1), wts would potentially look like
   * wts = [ F(i=1,n=-1), F(i=1,n=-2), ..., F(i=1,n=1), F(i=1,n=2), ... ]
   * except not really, because we only count not-small entries (>1e-8)
   *
   * For the second oscillator (i=2) we add in cross terms wherever is nice
   * wts = [ F(i=1,n=-1)*F(i=2,n=-1), F(i=1,n=-1)*F(i=2,i=-2), ... ]
   * given, of course, that they're reasonably large enough
   *
   * So hopefully after calling this enough, wts will be quite populated with
   * a lot of F(i_a,n_b)*F(i_c,n_d)*F(i_e,n_f)*.... terms. 
   */

  /* There are 50 entries in bplus and bminus (not all of them necessarily
   * nonzero). Loop through them. 
   */
  for ( auto k = 0; k < 50; ++k ){
    // It doesn't really make sense for b_minus_or_plus[k] to be negative
    // right? That's definitely wrong, right? It can be 0 just not < 0
    if ( b_minus_or_plus[k] <= 0 ){ return; } 
    for ( auto m = 0; m < nn; ++m ){
      if ( wtn[m] * b_minus_or_plus[k] >= 1e-8 and n < int(bes.size())-1 ){
        n += 1;
        bes[n] = ben[m] + pos_or_neg * (k+1) * beta_i;
        wts[n] = wtn[m] * b_minus_or_plus[k];
      }
    }
  }
}


auto oscillatorLoop( const double& alpha, std::vector<double>& debyeWaller, 
  std::vector<double>& ar, std::vector<double>& wts, std::vector<double>& bes, 
  std::vector<double>& betaVals, double& tbart, std::vector<double>& 
  t_eff_consts, const double& temp ){
  /* alpha          --> yup
   * debyeWaller       --> weight / ( tanh( 0.5 * energy / tev ) * energy / tev )
   *                    --defined in Eq. 538, evaluated in prepareParams.h
   * ar             --> weight / ( sinh( 0.5 * energy / tev ) * energy / tev )
   *                    --bessel arg from Eq. 537, evaluated in prepareParams.h
   *                    --input into discre, will be sent in from leapr
   * wts            --> blank vector with maxdd = 500 entries
   * bes            --> blank vector with maxdd = 500 entries
   * betaVals       --> energy / tev
   * a              --> alpha index
   * wt             --> tbeta
   * tbart          --> T_eff / temp  ( t_eff_vec[itemp] / temp_vec[itemp] )
   * oscEnergiesWeights --> vector of delta function energies/weights in tuples
   * t_eff_consts   --> 0.5 * weight * energy / tanh( 0.5 * energy / tev )
   *                    --defined in Eq. 544, evaluated in prepareParams.h
   *                    --in leapr.f90, this is called dist
   */

  std::vector<double> ben(wts.size(), 0), wtn(wts.size(), 0); wtn[0] = 1.0;

  double x, bzero;
  int n = 0, nn = 0;

  // Loop over all oscillators
  for ( auto i = 0; i < int(betaVals.size()); ++i ){
    nn = n + 1;
    x              = alpha*ar[i];
    //             = alpha*weight / (sinh(0.5*energy/tev) * energy/tev)
    //             = argument for bessel function in Eq. 537

    /* bfact populates bplus and bminus with A_in terms from Eq. 537.
     * The nth entry of bplus or bminus corresponds to a specific alpha and 
     * i value. The bzero output is either 
     *                 I0(x)*e^(-alpha*debyeWaller)
     * or 
     *                I0(x)*e^(-alpha*debyeWaller+x)
     * depending on the size of x.
     */
    std::vector<double> bminus (50,0), bplus(50,0);
    bzero = bfact( x, alpha*debyeWaller[i], betaVals[i], bplus, bminus );
    
    // do convolution for delta function
    n = 0;
    for ( auto m = 0; m < nn; ++m ){
      if ( (ben[m] <= 0 or wtn[m]*bzero >= 1e-8) and n < int(wts.size())){
        bes[n] = ben[m];
        wts[n] = wtn[m]*bzero;
        // Why are we multiplying by bzero? Because that's the A_(i,n) term for
        // n = 0, which we need to sum over. But notice that in bfact we 
        // populate a bplus and a bminus, and are kind of implicitly leaving
        // n=0 out. So we're accounting for this here. In the manual (pg. 713)
        // it states that we want to worry about the n=0 term first, so that's
        // what we're kind of doing. 
        n += 1;
      }
    }
    n -= 1;

    // Read the description for posNegTerms to get a better feel for this.
    // Basically we're going to be populating wts with A_i,n terms muliplied
    // by each other, for many different i and n. 
    posNegTerms(n, betaVals[i], bminus, wts, wtn, bes, ben, nn, -1);// negative
    posNegTerms(n, betaVals[i], bplus,  wts, wtn, bes, ben, nn, 1 );// positive

    // Copy first n entries of permanent array into our temporary arrays
    for ( auto m = 0; m <= n; ++m ){
      ben[m] = bes[m];
      wtn[m] = wts[m];
    }

    // Effective temperature is amended, this ( kind of ) follows Eq. 544.
    tbart += t_eff_consts[i] / ( kb * temp );


  }   
  return n; 

}


