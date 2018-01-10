#include <iostream>
#include <vector>
#include "coher_util/formf.h"
#include "coher_util/hexLatticeFactors.h"
#include "coher_util/bccLatticeFactors.h"
#include "coher_util/fccLatticeFactors.h"
#include "coher_util/end.h"





auto coher( int iel, int npr, int maxb, std::vector<double>& b, 
  double& emax ){

  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   *
   * inputs :iel tells you the material composition, npr tell you the number of
   * primary atoms, the b vector is the bragg edges vector, nbe is the number 
   * of edges, maxb and emax are also there
   */
  int i,j,k,imax,ifl,nw,nbe;
  double amne,econ,tsqx,a,c,amsc,scoh,c1,c2,recon,scon,wint,t2,
    ulim,phi,w1,w2,w3,tsq,tau,w,f,x,bel,be,bs,
    gr1 = 2.4573e-8,
    gr2 = 6.700e-8,
    gr3 = 12.011,
    gr4 = 5.50,
    be1 = 2.2856e-8,
    be2 = 3.5832e-8,
    be3 = 9.01,
    be4 = 7.53,
    beo1 = 2.695e-8,
    beo2 = 4.39e-8,
    beo3 = 12.5,
    beo4 = 1.0,
    al1 = 4.04e-8,
    al3 = 26.7495,
    al4 = 1.495,
    pb1 = 4.94e-8,
    pb3 = 207.,
    pb4 = 1.,
    fe1 = 2.86e-8,
    fe3 = 55.454,
    fe4 = 12.9,
    toler = 1.e-6,
    eps = .05,
    amassn = 1.008664904,
    amu = 1.6605402e-24,
    ev = 1.60217733e-12,
    hbar = 1.05457266e-27;

  // initialize.
  amne = amassn * amu;          // Mass of neutron in g ( 1.6749286E-24 )
  econ = ev*8*(amne/hbar)/hbar;
  recon = 1/econ;
  tsqx = econ/20;

  // For hexagonal materials the lattice is described by two constants, a and c,
  // which are defined below. Among other things this is used to define the 
  // reciprocal lattice vector lengths, which is done using tausq. The other
  // reciprocal lattice vector lengths require only one constant.
  if (iel == 1){
     // graphite constants.
     a=gr1;
     c=gr2;
     amsc=gr3;
     scoh=gr4/npr;
  }
  else if (iel == 2) {
     //  beryllium constants
     a=be1;
     c=be2;
     amsc=be3;
     scoh=be4/npr;
  }
  else if (iel == 3) {
     //  beryllium oxide constants
     a=beo1;
     c=beo2;
     amsc=beo3;
     scoh=beo4/npr;
  }
  else if (iel == 4) {
     // aluminum constants
     a=al1;
     amsc=al3;
     scoh=al4/npr;
  }
  else if (iel == 5) {
     // lead constants
     a=pb1;
     amsc=pb3;
     scoh=pb4/npr;
  }
  else if (iel == 6) {
     // iron constants
     a=fe1;
     amsc=fe3;
     scoh=fe4/npr;
  }

  scon = scoh*16*M_PI*M_PI;
  if (iel < 4) {
     c1=4/(3*a*a);
     c2=1/(c*c);
     scon /= (2*a*a*c*sqrt(3)*econ);
  }
  else if (iel == 4 or iel == 5) {
     c1=3/(a*a);
     scon/=(16*a*a*a*econ);
  }
  else if (iel == 6) {
     c1=2/(a*a);
     scon/=(8*a*a*a*econ);
  }
  wint=0;
  t2=hbar/(2*amu*amsc);
  ulim=econ*emax;
  ifl=1;
  nw=maxb;


  if ( iel < 4 ){
    // compute lattice factors for hexagonal lattices
    int i1m = a*sqrt(phi) + 1;
    imax = hexLatticeFactors( a, tsq, c1, c2, iel, nw, tsqx, b, ifl,  
    i, wint, t2, ulim, imax, c, i1m );
    k = imax + 1;
  }

  if ( iel < 6 ){
    // compute lattice factors for fcc lattices
    imax = fccLatticeFactors( iel, b, ifl, w, nw, t2, c1, wint, ulim, a );
    k = imax + 1;
  } 

  if ( iel == 6 ){
    // compute lattice factors for bcc lattices
    imax = bccLatticeFactors( ulim, b, ifl, wint, t2, iel, a, c1 );
    k = imax + 1;
  }

  
  end( ifl, b, k, recon, maxb, toler, scon, nw, ulim, imax );
 
  nbe = 1;
  return nbe;

}

