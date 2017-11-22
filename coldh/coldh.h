#include <iostream>
#include <vector>


auto coldh( int itemp, double temp, double tev, double sc, int ncold,
    double trans_weight, double tbeta, const std::vector<double>& tempf,
    const std::vector<double>& tempr, double scaling, 
    std::vector<double>& alpha ){
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

  
  int law = ncold + 1;
  double de, x, amassm, bp, sampc, sampi, wt, tbart;



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
  }
    


















   
}   
