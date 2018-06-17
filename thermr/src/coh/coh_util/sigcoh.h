#include "coh/coh_util/sigcoh_util/form.h"
#include "coh/coh_util/sigcoh_util/terp.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include <iostream>

bool finish( int& k, std::vector<double>& wrk, const double& f, int nw,
  const double& tau_sq ){
  k = k + 1;
  if ((2*k) > nw) { 
    std::cout << "storage exceeded" << std::endl;
    return true;
  }
  wrk[2*k-2] = tau_sq; wrk[2*k-1] = f;
  return false;
}


auto tausq( int m1, int m2, int m3, double c1, double c2 ){
  /* This function computes the value tau^2, which is defined in the General
   * Atomics HEXSCAT documentation. On pg. 62 of the HEXSCAT pdf document
   * (which is in HEXSCAT Appendix, Section 1: Formulation), we have that 
   *
   *     tau^2/4pi^2  = (4/3a^2) * ( l1^2 + l2^2 + l1*l2 ) + l3 / c^2
   *
   */
  return (c1*(m1*m1+m2*m2+m1*m2)+c2*(m3*m3))*4*M_PI*M_PI;
}

void swapVals( double& a, double& b ){
  double c = a; a = b; b = c;
}
  
auto computeCrossSections( double e, std::vector<double>& fl, 
    std::vector<double>& s, double emax, double scon, double recon, int nl, 
    std::vector<double> p, int k ){
  // compute cross sections at this energy
   double elim;
   for ( int il = 0; il < nl; ++il ){
      s[il]=0;
   }
   int last=0;

   for ( int i = 1; i <= k; ++i ){
      double tau_sq=fl[2*i-2];
      elim=tau_sq*recon;
      // if (elim >= e) exit
      double f=fl[2*i-1];
      if (e > emax) f=0;
      double u=1-2*elim/e;
      // u here is equal to fl for l = 1 (P1 component).
      // This is defined in the General Atomics HEXSCAT paper, in Part 1 
      // Formulation. If l == 0, fl = 1. But if l == 1, then
      //            fl = 1 - tau^2 lambda^2 / 8 pi^2
      //    which simplifies to 
      //                 1 - tau^2 hbar^2 / 4 m_n E

      int lmax=nl-1;
      legndr(u,p,lmax);
      for ( int il = 0; il < nl; ++il ){
         s[il]=s[il]+f*p[il];
      }
      if (i == k) last=1;
   }
   for ( int il = 0; il < nl; ++il ){
      s[il]=s[il]*scon/e;
   }
   if (last == 1 or elim > emax ) { elim=emax; }
}




auto sigcoh( double e, double& enext, std::vector<double> s, int nl, int lat, 
  double temp, double emax, int natom, std::vector<double>& fl, 
  std::vector<double>& p, int k, double scon ){
 /*-------------------------------------------------------------------
  * Compute the first nl Legendre components of the coherent scatter-
  * ing at energy e from lattice type lat.  Here enext is the next
  * Bragg edge.  Initialize if e=0.  A list of reciprocal lattice
  * shells and weights is precomputed and stored for use at all e.
  * Long, closely-spaced shells are grouped together to speed up the
  * calculation.
  *       lat=1  graphite
  *       lat=2  be
  *       lat=3  beo
  *       lat=10 read from endf6
  * nl returns no. of Bragg edges on initialization call.
  *-------------------------------------------------------------------
  */
  int nw, i1m, l1, i2m, i3m;
  double amne, econ, tsqx, a, c, amsc, scoh, wal2, wint, w1, w2, w3, tau_sq, 
    tau, w, f, t2, ulim, phi, c1, c2;
  
  // These are lattice factors. Apparently they were borrowed directly from
  // HEXSCAT code. 
  double gr1 = 2.4573e-8, // http://www.phy.ohiou.edu/~asmith/NewATOMS/HOPG.pdf
         gr2 = 6.700e-8,  // http://www.phy.ohiou.edu/~asmith/NewATOMS/HOPG.pdf 
         be1 = 2.2856e-8, // http://periodictable.com/Properties/A/LatticeConstants.html
         be2 = 3.5832e-8, // http://periodictable.com/Properties/A/LatticeConstants.html 
         beo1 = 2.695e-8, // https://link.springer.com/chapter/10.1007%2F10681719_737
         beo2 = 4.39e-8;  // https://link.springer.com/chapter/10.1007%2F10681719_737
                          // II-VI and I-VII Compounds; Semimagnetic Compounds

  // These are masses
  double gr3 = 12.011e0, 
         be3 = 9.01e0, 
         beo3 = 12.5e0;  // Mass of BeO is actually 25, but apparently we 
                         // divide by 2 because I suppose avg mass per atom
                         
  // These are the characteristic coherent cross sections for hte material.
  // These first appear in Eq. 222 on pg. 166.
  double gr4 = 5.50, // pg. 18 Neutron Physics Karl-Heinrich Beckurts, Karl Wirtz
         be4 = 7.53, // pg. 18 Neutron Physics Karl-Heinrich Beckurts, Karl Wirtz
         beo4 = 1.0;

  double cw = 0.658173e-15, hbar = 1.05457266e-27, amu = 1.6605402e-24, 
         amassn = 1.008664904, ev = 1.60217733e-12;
  
  // eps is the current grouping factor, 5%. This is used to lump together 
  // multiple tau values, so as to save storage and run time.
  double eps = 0.05;

  // save k,recon,scon


 /* If this is the first entry (E=0) for an ENDF-III type material, the 
  * appropriate lattice constants are selected and the Debye-Waller coefficient 
  * is obtained for the desired temperature by interpolation. Then the 
  * reciprocal lattice wave vectors and structure factors are computed, 
  * sorted into shells, and stored for later use.
  */
  amne = amassn * amu;     // mass of neutron in grams
  econ = ev * 8 * ( amne / hbar ) / hbar;
  tsqx = econ / 20;


 /* If energy is greater than zero, the stored list is used to compute the 
  * cross section. For ENDF-6 format materials, the initialization step is 
  * used to organize the data already read from MF=7/MT=2 by rdelas, and 
  * subsequent entries are used to compute the cross section.
  */
  if (e > 0) {
    double recon = 1.0/econ;
    computeCrossSections( e, fl, s, emax, scon, recon, nl, p, k );
    return std::vector<double> {};
  }
 

  // Temperatures interpolated over when trying to get correct Debye-Waller
  // Coefficient.
  std::vector<double> tmp {296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000};

  // Debye Waller to interpolate over
  // This appears in Eq. 220 of the manual, as 'W' in the exponential term
  // This also appears in the General Atomics HEXSCAT report as the integral
  // shown below (as wal2 description)

  std::vector<double> dwf ( 10 );
  if (lat == 1){       // GRAPHITE
    a = gr1; c = gr2; amsc = gr3; scoh = gr4/natom; 
    dwf = { 2.1997, 2.7448, 3.2912, 3.8510, 4.4210, 4.9969, 6.1624, 7.3387, 
      9.6287, 11.992 };      
  }
  else if (lat == 2) { // BERYLLIUM
    a = be1; c = be2; amsc = be3; scoh = be4/natom;
    dwf = { 3.16663, 3.88842, 4.62944, 5.40517, 6.19880, 7.0042, 8.63665, 
      10.2865, 0, 0 };  
  }
  else if (lat == 3){  // BERYLLIUM OXIDE
    a = beo1; c = beo2; amsc = beo3; scoh = beo4/natom;
    dwf = {  2.153, 2.6374, 3.1348, 3.6513, 4.1798, 4.7164, 5.8052, 6.9068, 
      0, 0 };
  } 
  else {
    std::cout << "OH NO! Error over here. Illegal lat value" << std::endl;
    throw std::exception();
  }


  // wal2 is (supposed to be) equal to
  //
  //       1   |' f(w)        
  //      ---  |  ----  coth( w / kb T ) dw
  //       M  _|   w   
  //
  // this seems to be the case from pg. 5 of the General Atomics HEXSCAT code,
  // in the table of descriptions for input values into the original HEXSCAT.
  //
  // In this case, the main equation we're trying to compute ( Eq. 1 in the
  // GA HEXSCAT documentation, on pg. 1 ) can be computed using an exponential 
  // term 
  //
  //                exp[ ( -hbar^2 tau^2 / 2 ) * wal2 ]
  //

  wal2 = terp(tmp,dwf,temp,2);

  c1 = 4.0 / ( 3.0 * a * a );
  c2 = 1.0 / ( c * c );
  scon = scoh * ( 16.0 * M_PI*M_PI )/( 2.0 * a * a * c * sqrt(3) * econ );
  wint = cw * amsc * wal2;
  t2 = hbar / ( 2.0 * amu * amsc );

  // This is the tau^2 value that corresponds to the maximum considered energy.
  // Calculated according to Eq. 223.
  ulim = econ * emax;

  nw = 10000;
  std::vector<double> wrk(nw,0.0);


  // compute and sort lattice factors.
  k = 0;

  // phi = ( tau_max/2pi )^2
  phi = ulim / ( 4 * M_PI * M_PI ); 

  // l1 --> 0 : a * tau_max / 2pi + 1
  // l1max = alpha * sqrt(phi), on pg 63 of the HEXSCAT document pdf. 
  i1m = a * sqrt(phi) + 1;

  for ( int l1 = 0; l1 < i1m; ++l1 ){
    i2m = 0.5 * ( l1 + sqrt( 3 * ( a * a * phi - l1 * l1 ) ) ) + 1;

    for ( int l2 = l1; l2 < i2m; ++l2 ){
      i3m = c * sqrt(phi - c1*( l1*l1 + l2*l2 - l1*l2 )) + 1;

      for ( int l3 = 0; l3 < i3m; ++l3 ){

        // w1 is equal to M1 on pg 3 of General Atomics HEXSCAT appendix. 
        w1 = (l1 == l2) ? 1 : 2;

        // w2 is equal to M2 on pg 3 of General Atomics HEXSCAT appendix. 
        if      (l1 == 0 and l2 == 0) { w2 = 0.5; }  // First l1, l2 iteration
        else if (l2 == 0)             { w2 = 1;   }  // Any l1, first l2
        else                          { w2 = 2;   }  
        
        // w3 is equal to M3 on pg 3 of General Atomics HEXSCAT appendix. 
        w3 = (l3 == 0) ? 1 : 2;
         
        // We consider l2 and -l2 because of Eq. 4 on pg. 4 of the General 
        // Atomics HEXSCAT appendix.
        for ( double&& l2: { l2, -l2 } ){
          tau_sq = tausq(l1,l2,l3,c1,c2);

          if (tau_sq > 0 and tau_sq <= ulim ){
            tau = sqrt(tau_sq);

            // w1*w2*w3 --> weighting factor M, as dfined on pg. 3 of the 
            // General Atomics HEXSCAT appendix
            w = exp(-tau_sq*t2*wint)*w1*w2*w3/tau;
            f = w * form( lat, l1, l2, l3 );

            if (k > 0 and tau_sq > tsqx) {
              for ( int i = 1; i <= k; ++i ){

                // As tau_i gets large, the values of tau_i get more and more 
                // closely spaced together, so a range of tau values can be 
                // lumped together to give a single effective tau_i and f_i.
                // This uses a 5% (eps) grouping factor.
                if (tau_sq < wrk[2*i-2] or tau_sq >= (1+eps)*wrk[2*i-2]){
                  if ( i == k ){ 
                    if ( finish( k, wrk, f, nw, tau_sq ) ){ return wrk; }
                    break;
                  }
                }
                else {                        // because got rid of continue
                  wrk[2*i-1]=wrk[2*i-1]+f;    // statement
                  break;
                }
              }
            }
            else {
              if ( finish( k, wrk, f, nw, tau_sq ) ){ return wrk; }
            }
          } 
        }
      }
    } 
  }

  for ( int i = 1; i <= k - 1; ++i ){
    for ( int j = i + 1; j <= k; ++j ){
      if (wrk[2*j-1-1] < wrk[2*i-1-1]) {
        swapVals( wrk[2*j-2], wrk[2*i-2] );
        swapVals( wrk[2*j-1], wrk[2*i-1] );
      } 
    }
  }
  k += 1;
  wrk[2*k-2]=ulim;
  wrk[2*k-1]=wrk[2*k-2-1];
  nw=2*k;
  enext=wrk[1-1]/econ;
  nl=k;
  return wrk;
  // return wrk as fl

}

















