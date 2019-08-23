#include "generalTools/interpolate.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>

template <typename Float, typename Range>
auto skold( Float cfrac, Float tev,const Range& alpha, const Range& beta, 
  const Range& skappa, Float awr, Float dka, Range& sab ){

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
   * sab  : S(a,b) that we've processed so far. Should only have incoherent
   *           contributions at this point I think.
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * SymSab is the modified S(a,b), following the above description of 
   *   the Skold Approximation. 
   */

  int a2;
  Float sk, ap;
  Range scoh( alpha.size() );
  Range kappa2Grid = ranges::view::iota(0,int(skappa.size()));
  Range kappaGrid = kappa2Grid 
                  | ranges::view::transform([delta=dka](auto& x){return delta*x;});

  // wave number [inverse angstroms]
  Float tempJoules = tev*ev;
  Float wavenVariables = 1.0e-10 * sqrt(2*awr*massNeutron*tempJoules) / hbar;
  Range waven = alpha | ranges::view::transform([wavenVariables](Float a){
                           return wavenVariables*sqrt(a);});

  // apply the skold approximation
  for ( size_t b = 0; b < beta.size(); ++b ){    
    for ( size_t a = 0; a < alpha.size(); ++a ){

      // Interpolate to find the waven value in the skappa input 
      sk = interpolate( skappa, waven[a], kappaGrid, 1.0 );
      ap = alpha[a] / sk;
      for ( a2 = 0; a2 < int(alpha.size()); a2++ ){
        if (ap < alpha[a2]){ break; }
      }
      if (a2 == int(alpha.size())){ a2 -= 1; }
      if (a2 == 0) { a2 = 1; }

      // Interpolate to calculate S(a',b) where a' = a/S(k)
      scoh[a] = terp1( alpha[a2-1], sab[b+(a2-1)*beta.size()], alpha[a2], 
             sab[b+a2*beta.size()], ap, 5 ) * sk;
    } // alpha loop

    // Amend the existing scattering law by combining a piece of it with a 
    // piece of the coherent scattering interactions. 
    for ( size_t a = 0; a < alpha.size(); ++a ){
      sab[b+a*beta.size()] = (1-cfrac)*sab[b+a*beta.size()]+cfrac*scoh[a];
    } // alpha loop
  } // beta loop
}

