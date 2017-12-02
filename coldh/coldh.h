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
  double deh = 0.0147;
  double ded = 0.0074; 
  double amassh = 3.3465E-24;
  double amassd = 6.69E-24;
  double pmass = 1.6726231E-24;
  double dmass = 3.343568E-24;
  double hbar = 1.05457266e-27;
  double sampch = 0.356;
  double sampcd = 0.668;
  double sampih = 2.526;
  double sampid = 0.403;
  double therm = 0.0253;

  int nbx; 
  int law = ncold + 1;
  int maxbb = 2 * beta.size() + 1;
  double de, x, amassm, bp, sampc, sampi, wt, tbart;


  std::vector<double> exb(maxbb, 0.0 );
  std::vector<double> betan(nbeta, 0.0 );
  std::vector<double> bex(maxbb, 0.0 );
  std::vector<double> rdbex(maxbb, 0.0 );
 

  if ( law > 3 ){
    de = ded;
    amassm = amassd;
    bp = hbar * sqrt(2/(ded*eV*dmass)) / ( 2 * angst ); 
    sampc = sampcd;
    sampi = sampid;
  } else {
    de = deh;
    amassm = amassh;
    bp = hbar/2*sqrt(2/(deh*eV*pmass))/angst;
    sampc = sampch;
    sampi = sampih;
  }

  x = de / tev;
  wt = trans_weight + tbeta;
  tbart = tempf[itemp] / tempr[itemp];

  //std::cout << bp << std::endl;


  for ( auto a = 0; a < alpha.size(); ++a ){
    double al = alpha[a]*scaling;
    double alp = wt * al;
    double waven = angst * sqrt( amassm * tev * eV * al ) / hbar;
    double y = bp * waven;
    double sk = terpk( ska, dka, waven );

    // spin-correlation factors
    double swe, swo, snorm;
    //std::cout << law << std::endl;
    if (law == 2){ 
      swe=sampi*sampi/3;
      swo=sk*sampc*sampc+2*sampi*sampi/3;
    } else if ( law == 3 ){
      swe=sk*sampc*sampc;
      swo=sampi*sampi;
    } else if ( law == 4 ){
      swe=sk*sampc*sampc+5*sampi*sampi/8;
      swo=3*sampi*sampi/8;
    } else if (law == 5){ 
      swe=3*sampi*sampi/4;
      swo=sk*sampc*sampc+sampi*sampi/4;
    }
    snorm=sampi*sampi+sampc*sampc;
    swe=swe/snorm;
    swo=swo/snorm;

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

    betaLoop( betan, rdbex, bex, sex, al, wt, tbart, x, y, swe, swo, itemp, 
       nbx, a, law, sym_sab, sym_sab_2 );

   
  }
 
    return;
   
}   
