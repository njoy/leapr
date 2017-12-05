#include <iostream>
#include <vector>
#include "coldh_util/terpk.h"
#include "coldh_util/betaLoop.h"
#include "../discre/discre_util/bfill.h"
#include "../discre/discre_util/exts.h"


auto coldh( int itemp, double temp, double tev, double sc, int ncold,
    double trans_weight, double tbeta, const std::vector<double>& tempf,
    const std::vector<double>& tempr, double scaling, 
    const std::vector<double>& alpha, const std::vector<double>& beta, 
    double& dka, std::vector<double>& ska, int nbeta, int lat,
    std::vector<std::vector<std::vector<double>>>& sym_sab,
    std::vector<std::vector<std::vector<double>>>& sym_sab_2 ){
  /* Convolve current scattering law with discrete rotational modes for ortho
   * or para hydrogen / deuterium. The discrete modes are calculated using 
   * formulas of Young and Koppel for vibrational ground state with coding 
   * based on contributions from Robert (Grenoble) and Neef (Julich). The 
   * approach of using solid/diffusive modes with discrete rotations is based
   * on the work of Keinert and Sax. Note that the final S(a,b) is not 
   * symmetric in beta
   */

  double angst = 1.0e-8;
  double eV = 1.60217733e-12;
  double massH = 1.6726231E-24;
  double massD = 3.343568E-24;
  double hbar = 1.05457266e-27;
  double therm = 0.0253;

  int nbx; 
  int maxbb = 2 * beta.size() + 1;
  double de, x, amassm, bp, scatLenCoh, scatLenIncoh, wt, tbart;


  std::vector<double> exb(maxbb, 0.0 );
  std::vector<double> betan(nbeta, 0.0 );
  std::vector<double> bex(maxbb, 0.0 );
  std::vector<double> rdbex(maxbb, 0.0 );
 

  // Either Ortho Deuterium or Para Deuterium 
  if ( ncold > 2 ){
    de = 0.0074;
    amassm = 6.69E-24;
    bp = hbar * sqrt( 2 /( de*eV*massD ) ) / ( 2 * angst ); 
    scatLenCoh = 0.668;
    scatLenIncoh = 0.403;
  } 
  // Either Ortho Hydrogen or Para Hydrogen
  else {
    de = 0.0147;
    amassm = 3.3465E-24;
    bp = hbar * sqrt( 2 /( de*eV*massH ) ) / ( 2 * angst );
    scatLenCoh = 0.356;
    scatLenIncoh = 2.526;
  }

  x = de / tev;
  wt = trans_weight + tbeta;
  tbart = tempf[itemp] / tempr[itemp];



  for ( auto a = 0; a < alpha.size(); ++a ){
    double al = alpha[a]*scaling;
    double alp = wt * al;
    double waven = angst * sqrt( amassm * tev * eV * al ) / hbar;
    double y = bp * waven;
    double sk = terpk( ska, dka, waven );

    // spin-correlation factors
    double evenSum, oddSum;
    // -----------------------------------------------------------------------
    // Thie is meant to recreate the table on pg. 662 of the manual, where we
    // get the A (even) and B (odd) terms for the summation in Eq. 567.
    // -----------------------------------------------------------------------
    // Ortho Hydrogen
    if (ncold == 1){ 
      evenSum = scatLenIncoh * scatLenIncoh/3;
      oddSum = sk * scatLenCoh * scatLenCoh + 2 * scatLenIncoh * scatLenIncoh / 3;
    } 
    // Para Hydrogen
    else if ( ncold == 2 ){
      evenSum = sk * scatLenCoh * scatLenCoh;
      oddSum = scatLenIncoh * scatLenIncoh;
    } 
    // Ortho Deuterium
    else if ( ncold == 3 ){
      evenSum = sk * scatLenCoh * scatLenCoh + 5 * scatLenIncoh * scatLenIncoh / 8;
      oddSum = 3 * scatLenIncoh * scatLenIncoh / 8;
    } 
    // Para Deuterium
    else if ( ncold == 4){ 
      evenSum = 3 * scatLenIncoh * scatLenIncoh / 4;
      oddSum = sk * scatLenCoh * scatLenCoh + scatLenIncoh * scatLenIncoh / 4;
    }

    evenSum = evenSum / (scatLenIncoh * scatLenIncoh + scatLenCoh * scatLenCoh);
    oddSum = oddSum / (scatLenIncoh * scatLenIncoh + scatLenCoh * scatLenCoh);

    // prepare arrays for sint
    
   if (a == 0){ 
      for ( auto b = 0; b < nbeta; ++b ){
          double be=beta[b];
          if (lat == 1){ be = be * therm / tev; }
          exb[b] = exp(-be/2);
          betan[b] = be;
      } 
      nbx = bfill(bex,rdbex,betan);
    }
    std::vector<double> input ( beta.size(), 0.0 ); 
    for ( auto b = 0; b < beta.size(); ++b ){
      input[b] = sym_sab[a][b][itemp];
    }
    auto sex = exts( input, exb, betan );

    betaLoop( betan, rdbex, bex, sex, al, wt, tbart, x, y, evenSum, oddSum, itemp, 
       nbx, a, ncold, sym_sab, sym_sab_2 );

   
  }
 
    return;
   
}   
