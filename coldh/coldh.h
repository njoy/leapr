#include <iostream>
#include <vector>


auto coldh( int itemp, double temp, double tev, double sc, int ncold ){
  /* Convolve current scattering law with discrete rotational modes for ortho
   * or para hydrogen / deuterium. The discrete modes are calculated using 
   * formulas of Young and Koppel for vibrational ground state with coding 
   * based on contributions from Robert (Grenoble) and Neef (Julich). The 
   * approach of using solid/diffusive modes with discrete rotations is based
   * on the work of Keinert and Sax. Note that the final S(a,b) is not 
   * symmetric in beta
   */

   double angst = 1.0e-8;
   double eV = 1.60217733e-12
   double deh = 0.0147;
   double ded = 0.0074; 
   double amassh = 3.3465E-24;
   double amassd = 6.69E-24;
   double pmass = 1.6726231E-24;
   double dmass = 3.343568E-24;
  
   int law = ncold + 1;
   double de = deh;
   if ( law > 3 ){ de = ded; }
   double x = de / tev;
  
   double amassm = amassh;
   if ( law > 3 ){ amassm = amassd; }
   double bp = hbar/2*sqrt(2/(deh*ev*pmass))/angst;
   if ( law > 3 ){ amassm = amassd; }
   
}   
