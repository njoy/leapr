#include <iostream>
#include <vector>
#include "coldh_util/terpk.h"
#include "coldh_util/betaLoop.h"
#include <unsupported/Eigen/CXX11/Tensor>


auto coldh( int itemp, const double& temp, double tev, int ncold,
    double trans_weight, double tbeta, const std::vector<double>& tempf,
    double scaling, 
    const std::vector<double>& alpha, const std::vector<double>& beta, 
    double& dka, std::vector<double>& ska, int nbeta, int lat, bool free, 
    Eigen::Tensor<double,3>& sym_sab,
    Eigen::Tensor<double,3>& sym_sab_2 ){
  /* Convolve current scattering law with discrete rotational modes for ortho
   * or para hydrogen / deuterium. The discrete modes are calculated using 
   * formulas of Young and Koppel for vibrational ground state with coding 
   * based on contributions from Robert (Grenoble) and Neef (Julich). The 
   * approach of using solid/diffusive modes with discrete rotations is based
   * on the work of Keinert and Sax. Note that the final S(a,b) is not 
   * symmetric in beta
   */

  double angst  = 1.0e-10,
         eV     = 1.60217733e-19,  // elementary charge [J] 
         mass_H = 1.6726231E-27,   // Mass of H in grams
         mass_D = 3.343568E-27,    // Mass of D in grams
         hbar   = 1.05457266e-34,  // planck constant [J*s]
         de, x, mass_H2_D2, bp, scatLenC, scatLenI, wt, tbart, therm = 0.0253;
  int nbx, maxbb = 2 * beta.size() + 1;


  std::vector<double> exb(maxbb, 0.0), betan(nbeta, 0.0), bex(maxbb, 0.0), 
    rdbex(maxbb, 0.0);
 

  // Either Ortho Deuterium or Para Deuterium 
  if ( ncold > 2 ){
    de = 0.0074;
    mass_H2_D2 = 6.69E-27;    // Mass of D2 in kg
    bp = hbar / ( angst * sqrt( 2*de*eV*mass_D ) ); 
    scatLenC = 0.668; // If you want to see where these values probably come
    scatLenI = 0.403; // from, consider looking at 
    // https://www.ncnr.nist.gov/resources/n-lengths/elements/h.html
  } 
  // Either Ortho Hydrogen or Para Hydrogen
  else {
    de = 0.0147;
    mass_H2_D2 = 3.3465E-27;  // Mass of H2 in kg
    // It seems like this is trying to be a/2, as defined on pg. 662, but when
    // I plug in values it seems to be 2 orders of magnitude off.
    bp = hbar / ( angst * sqrt( 2*de*eV*mass_H ) );
    scatLenC = 0.356;
    scatLenI = 2.526;
  }

  x = de / tev;
  wt = trans_weight + tbeta;   // Translational weight
  tbart = tempf[itemp] / temp; // Effective temperature



  for ( size_t a = 0; a < alpha.size(); ++a ){
    double al = alpha[a]*scaling;
    double waven = angst * sqrt( mass_H2_D2 * tev * eV * al ) / hbar;
    double y = bp * waven;

    // We interpolate S(kappa) to get the corresponding value that we will use
    // for the Vineyard approximation. 
    // We will replace a_c^2 with S(kappa)*a_c^2, as is stated on pg. 663.
    double sk = terpk( ska, dka, waven );

    // spin-correlation factors
    double evenSumConst, oddSumConst;
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

    
    // Ortho Hydrogen
    if (ncold == 1){ 
      evenSumConst =      scatLenI * scatLenI / 3;
      oddSumConst  = sk * scatLenC * scatLenC + 2 * scatLenI * scatLenI / 3;
    } 
    // Para Hydrogen
    else if ( ncold == 2 ){
      evenSumConst = sk * scatLenC * scatLenC;
      oddSumConst  =      scatLenI * scatLenI;
    } 
    // Ortho Deuterium
    else if ( ncold == 3 ){
      evenSumConst = sk * scatLenC * scatLenC + 5 * scatLenI * scatLenI / 8;
      oddSumConst  = 3  * scatLenI * scatLenI / 8;
    } 
    // Para Deuterium
    else {// if ( ncold == 4){ 
      evenSumConst = 3  * scatLenI * scatLenI / 4;
      oddSumConst  = sk * scatLenC * scatLenC + scatLenI * scatLenI / 4;
    }

    // Both the ortho and para scattering laws (Eq. 567-568) are multiplied by
    // a (4*pi/sigma_b) term (where sigma_b is the characteristic bound cross
    // section, 
    //    sigma_b = 4*pi*(Incoh Scattering Len^2 + Coh Scattering Len^2 ).
    // So to account for the 4*pi/sigma_b term, we divide by 
    //          Incoh Scattering Len^2 + Coh Scattering Len^2
    evenSumConst /= scatLenI * scatLenI + scatLenC * scatLenC;
    oddSumConst  /= scatLenI * scatLenI + scatLenC * scatLenC;

    // prepare arrays for sint
    
   if (a == 0){ 
      for ( int b = 0; b < nbeta; ++b ){
          double be=beta[b];
          if (lat == 1){ be = be * therm / tev; }
          exb[b] = exp(-be);
          betan[b] = be;
      } 
      nbx = bfill(bex,rdbex,betan);
    }
    std::vector<double> input ( beta.size(), 0.0 ); 
    for ( size_t b = 0; b < beta.size(); ++b ){
      input[b] = sym_sab(a,b,itemp);
    }
    auto sex = exts( input, exb, betan );

    betaLoop( betan, rdbex, bex, sex, al, wt, tbart, x, y, evenSumConst, 
      oddSumConst, itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );

   
  }
 
}   
