#include <iostream>
#include <vector>
#include "../coldh/coldh_util/terpk.h"
#include "skold_util/terp1.h"

auto skold( double cfrac, int itemp, double tev,
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<double>& skappa, double awr, 
  double dka, double scaling,
  std::vector<std::vector<std::vector<double>>>& symSab ){
  /* Overview 
   * ------------------------------------------------------------------------
   * The purpose of this is to apply the Skold approximation to add in the 
   * effects of intermolecular coherence. This results in a scattering law of
   *
   *               S(a,b) = (1-c)*S_inc(a,b) + (c)*S_coh(a,b)
   *
   * The Skold approximation to the coherent component is defined to be 
   *
   *              S_coh(a,b) = S( alpha/S(k0), b ) * S(k0)
   *
   * where S(kappa) is the static structure factor, and k0 is a kappa value in
   * units of inverse Angstroms. Check out the INDC Thermal Neutron Scattering
   * Data for the Moderator Materials H20 D20 .... by Mattes and Keinert April
   * 2005 for some information, or also could just got to Skold's original 
   * paper which is Small Energy Transfer Scattering of Cold Neutrons from
   * Liquid Argon, Eq. 2. Both also introduce Vineyard approximation.
   *
   *
   * Inputs
   * ------------------------------------------------------------------------
   * cfrac   : fraction that you use to average the incoherent and coherent
   *           contributions to the scattering law. 
   * itemp   : temperature index
   * tev     : temperature in eV
   * alpha   : alpha vector
   * beta    : beta vector
   * skappa  : S(kappa) vector, which is user input value from Card18 and 
   *           represents the static structure factor.
   * awr     : atomic weight ratio, calculated in leapr.cpp
   * dka     : spacing for the S(kappa) vector
   * scaling : value at which the alpha values should be scaled
   * symSab  : S(a,b) that we've processed so far. Should only have incoherent
   *           contributions at this point I think.
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * SymSab is the modified S(a,b), following the above description of 
   *   the Skold Approximation.
   */


  /* use skold approximation to add in the effects
   * of intermolecular coherence.
   */
  int kk;
  double sk, ap, waven, amassn = 1.008664904, 
         amu = 1.6605402e-24, bk = 8.617385e-5, hbar = 1.05457266e-27, 
         ev  = 1.60217733e-12;

  std::vector<double> scoh ( 1000, 0.0 );
  // apply the skold approximation
  for ( auto b = 0; b < beta.size(); ++b ){
    for ( auto a = 0; a < alpha.size(); ++a ){
      // Getting a value in units of inverse angstroms so that we can happily
      // interpolate it in our skappa table
      waven = 1.0e-8 * sqrt(2*awr*amassn*amu*tev*ev*alpha[a]*scaling)/hbar;
      // Interpolate to find the waven value in the skappa input 
      sk = terpk(skappa,dka,waven);
      ap = alpha[a] / sk;
      for ( auto a2 = 0; a2 < alpha.size(); ++a2 ){
        kk = a2;
        if (ap < alpha[a2]){ break; }
      }
      if (kk == 0) kk = 1;

      terp1( alpha[kk-1], symSab[kk-1][b][itemp], alpha[kk], 
             symSab[kk][b][itemp], ap, scoh[a], 5 );
      scoh[a] *= sk;
    }
    // Amend the existing scattering law by combining a piece of it with a 
    // piece of the coherent scattering interactions. 
    for ( auto a = 0; a < alpha.size(); ++a ){
      symSab[a][b][itemp] = (1-cfrac)*symSab[a][b][itemp]+cfrac*scoh[a];
    }
  }
}

