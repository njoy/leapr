#include <iostream>
#include <vector>
#include "coher_util/formf.h"
#include "coher_util/hexLatticeFactors.h"
#include "coher_util/bccLatticeFactors.h"
#include "coher_util/fccLatticeFactors.h"

auto sortLatticeFactors( int jmin, std::vector<double>& b, int ifl, int imax,
  int k, double& st, double& sf ){
  // sort lattice factors
  
  for ( auto i = 1; i < imax; ++i ){
    jmin=i+1;
    for ( auto j = jmin; j < k; ++j ){
      if (b[ifl+2*j-2-1] < b[ifl+2*i-2-1]) {
        st = b[ifl+2*i-2-1];
        sf = b[ifl+2*i-1-1];
        b[ifl+2*i-2-1] = b[ifl+2*j-2-1];
        b[ifl+2*i-1-1] = b[ifl+2*j-1-1];
        b[ifl+2*j-2-1] = st;
        b[ifl+2*j-1-1] = sf;
      }
    }
  }
}





auto coher( int lat, int natom, int nbe, int maxb, std::vector<double> b, 
  double emax ){
  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   */
  int i,j,k,imax,jmin,idone,ifl,nw,i1,i2,i3,l1,l2,l3,i2m,i3m;
  double time,twopis,amne,econ,tsqx,a,c,amsc,scoh,c1,c2,recon,scon,wint,t2,
    ulim,phi,w1,w2,w3,tsq,tau,w,f,x,st,sf,bel,be,bs;
  double gr1 = 2.4573e-8,
  gr2 = 6.700e-8,
  gr3 = 12.011e0,
  gr4 = 5.50e0,
  be1 = 2.2856e-8,
  be2 = 3.5832e-8,
  be3 = 9.01e0,
  be4 = 7.53e0,
  beo1 = 2.695e-8,
  beo2 = 4.39e-8,
  beo3 = 12.5e0,
  beo4 = 1.0e0,
  al1 = 4.04e-8,
  al3 = 26.7495e0,
  al4 = 1.495e0,
  pb1 = 4.94e-8,
  pb3 = 207.e0,
  pb4 = 1.e0,
  fe1 = 2.86e-8,
  fe3 = 55.454e0,
  fe4 = 12.9e0,
  twothd = 0.666666666667e0,
  sqrt3 = 1.732050808e0,
  toler = 1.e-6,
  eps = .05e0,
  pi = 3.14159265358979,
  amassn = 1.008664904,
  amu = 1.6605402e-24,
  ev = 1.60217733e-12,
  hbar = 1.05457266e-27;

  // initialize.
  twopis = (2*pi)*(2*pi);
  amne = amassn * amu;          // Mass of neutron in g ( 1.6749286E-24 )
  econ = ev*8*(amne/hbar)/hbar;
  recon = 1/econ;
  tsqx = econ/20;

  // For hexagonal materials the lattice is described by two constants, a and c,
  // which are defined below. Among other things this is used to define the 
  // reciprocal lattice vector lengths, which is done using tausq. The other
  // reciprocal lattice vector lengths require only one constant.
  if (lat == 1){
     // graphite constants.
     a=gr1;
     c=gr2;
     amsc=gr3;
     scoh=gr4/natom;
  }
  else if (lat == 2) {
     //  beryllium constants
     a=be1;
     c=be2;
     amsc=be3;
     scoh=be4/natom;
  }
  else if (lat == 3) {
     //  beryllium oxide constants
     a=beo1;
     c=beo2;
     amsc=beo3;
     scoh=beo4/natom;
  }
  else if (lat == 4) {
     // aluminum constants
     a=al1;
     amsc=al3;
     scoh=al4/natom;
  }
  else if (lat == 5) {
     // lead constants
     a=pb1;
     amsc=pb3;
     scoh=pb4/natom;
  }
  else if (lat == 6) {
     // iron constants
     a=fe1;
     amsc=fe3;
     scoh=fe4/natom;
  }

  scon = std::pow( scoh*4*M_PI, 2 );
  if (lat < 4) {
     c1=4/(3*a*a);
     c2=1/(c*c);
     scon /= (2*a*a*c*sqrt3*econ);
  }
  else if (lat == 4 or lat == 5) {
     c1=3/(a*a);
     scon/=(16*a*a*a*econ);
  }
  else if (lat == 6) {
     c1=2/(a*a);
     scon/=(8*a*a*a*econ);
  }
  wint=0;
  t2=hbar/(2*amu*amsc);
  ulim=econ*emax;
  ifl=1;
  nw=maxb;


  // compute lattice factors for hexagonal lattices
  int i1m = a*sqrt(phi) + 1;
  hexLatticeFactors( a, tsq, c1, c2, lat, nw, tsqx, b, ifl,  
  i, wint, t2, ulim, imax, c, i1m );

  // compute lattice factors for fcc lattices
  fccLatticeFactors( lat, b, ifl, w, nw, t2, c1, wint, ulim, a );

  // compute lattice factors for bcc lattices
  bccLatticeFactors( ulim, b, ifl, wint, t2, lat, a, c1 );

  // sort lattice factors
  sortLatticeFactors( jmin, b, ifl, imax, k, st, sf );
 

}

