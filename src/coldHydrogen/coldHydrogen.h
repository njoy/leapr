#include <iostream>
#include <vector>
#include "generalTools/interpolate.h"
#include "coldHydrogen_util/betaLoop.h"
#include <range/v3/all.hpp>
#include "generalTools/constants.h"
#include <tuple>


template <typename Float>
auto getSumConstants( int ncold, Float& scatLenI, Float& scatLenC, Float& sk ){
  using std::make_tuple; using std::pow;
  // Ortho Hydrogen
  if (ncold == 1){ 
    return make_tuple(                     pow(scatLenI,2)*0.333333,   // EVEN
                      pow(scatLenC,2)*sk + pow(scatLenI,2)*0.666667); }// ODD 

  // Para Hydrogen
  else if ( ncold == 2 ){
    return make_tuple(sk*pow(scatLenC,2),                              // EVEN
                         pow(scatLenI,2) ); }                          // ODD

  // Ortho Deuterium
  else if ( ncold == 3 ){
    return make_tuple(pow(scatLenC,2)*sk + pow(scatLenI,2)*0.625,     // EVEN
                                           pow(scatLenI,2)*0.375); }  // ODD

  // Para Deuterium
  else {   // ncold == 4{ 
    return make_tuple(                     pow(scatLenI,2)*0.75,      // EVEN
                      pow(scatLenC,2)*sk + pow(scatLenI,2)*0.25); }   // ODD
}




template <typename Float, typename Range>
auto coldHydrogen( Float tev, int ncold, Float transContinWeight, 
  const Range& alpha, const Range& beta, Float& dka, Range& ska, 
  bool free, Range& sab_1, Range& sab_2, const Float& tbart ){
  /* Convolve current scattering law with discrete rotational modes for ortho
   * or para hydrogen / deuterium. The discrete modes are calculated using 
   * formulas of Young and Koppel for vibrational ground state with coding 
   * based on contributions from Robert (Grenoble) and Neef (Julich). The 
   * approach of using solid/diffusive modes with discrete rotations is based
   * on the work of Keinert and Sax. Note that the final S(a,b) is not 
   * symmetric in beta
   */

  Float de, x, massMolecule, bp, scatLenC, scatLenI, wt;
  int nbx, maxbb = 2 * beta.size() + 1;

  Range exb(beta.size(), 0.0), bex(maxbb, 0.0), rdbex(maxbb, 0.0);

  // Either Ortho Deuterium or Para Deuterium 
  if ( ncold > 2 ){
    de = 0.0074;
    massMolecule = massD2;
    bp = hbar / ( 1e-10 * sqrt( 2*de*ev*massDeuterium ) ); 
    scatLenC = 0.668; // If you want to see where these values probably come
    scatLenI = 0.403; // from, consider looking at 
    // https://www.ncnr.nist.gov/resources/n-lengths/elements/h.html
  } 
  // Either Ortho Hydrogen or Para Hydrogen
  else {
    de = 0.0147;
    massMolecule = massH2;  
    // It seems like this is trying to be a/2, as defined on pg. 662, but when
    // I plug in values it seems to be 2 orders of magnitude off.
    bp = hbar / ( 1e-10 * sqrt( 2*de*ev*massHydrogen ) );
    scatLenC = 0.356;
    scatLenI = 2.526;
  }

  x = de / tev;
  //wt = trans_weight + tbeta;  
  wt = transContinWeight;

  auto xVals = ranges::view::iota(0,int(ska.size()))
             | ranges::view::transform([delta=dka](auto x){return Float(delta*x);}); 
  auto xyZipped = ranges::view::zip(xVals,ska);

  for ( int b = 0; b < int(beta.size()); ++b ){
      exb[b] = exp(-beta[b]*0.5);
  } 
  auto output = bfill(rdbex,beta);

  nbx = std::get<0>(output);
  for ( size_t i = 0; i < bex.size(); ++i ){
    bex[i] = std::get<1>(output)[i];
  }


  for ( size_t a = 0; a < alpha.size(); ++a ){
    Float waven = 1e-10 * sqrt( massMolecule * tev * ev * alpha[a] ) / hbar;
    Float y = bp * waven;

    // We interpolate S(kappa) to get the corresponding value that we will use
    // for the Vineyard approximation. 
    // We will replace a_c^2 with S(kappa)*a_c^2, as is stated on pg. 663.
    Float sk = interpolate( xyZipped, waven, 1.0 );

    // spin-correlation factors
    auto output = getSumConstants( ncold, scatLenI, scatLenC, sk );
    Float evenSumConst = std::get<0>(output), 
          oddSumConst  = std::get<1>(output);
    // -----------------------------------------------------------------------
    // This is meant to recreate the table on pg. 662 of the manual, where we
    // get the A (even) and B (odd) terms for the summation in Eq. 567.
    // -----------------------------------------------------------------------
    
    // Note that this differs slightly from Eq. 567-568, in that some of the 
    // terms are being multiplied by this sk term. This is S(kappa), which
    // is used to account for intermolecular coherence (if there is a 
    // correlation between positions of nearby molecules). S(kappa) is also 
    // called the static structure factor
    //
    // It seems that this is implicitly applying the Vineyard approximation,
    // which is detailed in Eq. 571 on pg. 663 of the NJOY manual. This 
    // seems to make sense, since if we wanted to account for intermolecular
    // coherence, we would either use Vineyard or Skold. Skold only gets 
    // invoked if we ask for it, AND if we don't use coldh. So the fact that
    // we're using this (and applying Vineyard) means that we are not able to
    // use Skold later. Which is good, because we wouldn't want to do it twice.


    // Both the ortho and para scattering laws (Eq. 567-568) are multiplied by
    // a (4*pi/sigma_b) term (where sigma_b is the characteristic bound cross
    // section, 
    //    sigma_b = 4*pi*(Incoh Scattering Len^2 + Coh Scattering Len^2 ).
    // So to account for the 4*pi/sigma_b term, we divide by 
    //          Incoh Scattering Len^2 + Coh Scattering Len^2
    evenSumConst /= scatLenI * scatLenI + scatLenC * scatLenC;
    oddSumConst  /= scatLenI * scatLenI + scatLenC * scatLenC;

    Range input ( beta.size(), 0.0 ); 
    for ( size_t b = 0; b < beta.size(); ++b ){
      input[b] = sab_1[b+a*beta.size()];
    }
    auto sex = extsCOLDH( input, exb, beta );

    betaLoop( beta, rdbex, bex, sex, alpha[a]*transContinWeight, tbart, x, y, 
      evenSumConst, oddSumConst, nbx, a, ncold, free, sab_1, sab_2 );

  }
 
}
    
