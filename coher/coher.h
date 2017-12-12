#include <iostream>
#include <vector>

auto coher( int lat, int natom, int nbe, int maxb, std::vector<double> b, 
  double emax ){
  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   */
  int i,j,k,imax,jmin,idone,ifl,i1m,nw,i1,i2,i3,l1,l2,l3,i2m,i3m;
  double time,twopis,amne,econ,tsqx,a,c,amsc,scoh,c1,c2,recon,scon,wint,t2,
    ulim,phi,w1,w2,w3,tsq,tau,w,f,x,st,sf,bel,be,bs
  gr1 = 2.4573e-8;
  gr2 = 6.700e-8;
  gr3 = 12.011e0;
  gr4 = 5.50e0;
  be1 = 2.2856e-8;
  be2 = 3.5832e-8;
  be3 = 9.01e0;
  be4 = 7.53e0;
  beo1 = 2.695e-8;
  beo2 = 4.39e-8;
  beo3 = 12.5e0;
  beo4 = 1.0e0;
  al1 = 4.04e-8;
  al3 = 26.7495e0;
  al4 = 1.495e0;
  pb1 = 4.94e-8;
  pb3 = 207.e0;
  pb4 = 1.e0;
  fe1 = 2.86e-8;
  fe3 = 55.454e0;
  fe4 = 12.9e0;
  twothd = 0.666666666667e0;
  sqrt3 = 1.732050808e0;
  toler = 1.e-6;
  eps = .05e0;
  zero = 0;

  // initialize.
  twopis = (2*pi)*(2*pi);
  amne = amassn * amu;
  econ=ev*8*(amne/hbar)/hbar;
  recon=1/econ;
  tsqx=econ/20;
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
  else {
    std::cout << "ERROR" << std::endl;
  }
}
