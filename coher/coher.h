#include <iostream>
#include <vector>
#include "coher_util/formf.h"
#include "coher_util/hexLatticeFactors.h"
#include "coher_util/bccLatticeFactors.h"

auto sqrtLatticeFactors( int jmin, std::vector<double>& b, int ifl, int imax,
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




auto end( int ifl, std::vector<double>& b, int k, int jmin, double recon, 
  double& bel, int nbe, int maxb, int j, double& be, double toler, double bs,
  double scon, double sf, double st, int nw, double ulim, int imax ){

  // Sort lattice factors
  // 220 continue
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
  k = k + 1;
  b[ifl+2*k-2-1] = ulim;
  b[ifl+2*k-1-1] = b[ifl+2*k-3-1];
  nw = 2 * k;

  // convert to practical units and combine duplicate bragg edges.
  bel=-1;
  j=0;
  for ( auto i = 1; i < k; ++i ){
    be = b[ifl+2*i-2-1]*recon;
    bs = b[ifl+2*i-1-1]*scon;
    if (be-bel < toler) {
      b[ifl+2*j-1-1] = b[ifl+2*j-1-1] + bs;
    }
    else {
      j = j+1;
      b[ifl+2*j-2-1] = be;
      b[ifl+2*j-1-1] = bs;
      bel=be;
    }
  }
  nbe = j;
  maxb = 2 * nbe;
}


auto coher( int lat, int natom, int nbe, int maxb, std::vector<double> b, 
  double emax ){
  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   */
  int i,j,k,imax,jmin,idone,ifl,i1m,nw,i1,i2,i3,l1,l2,l3,i2m,i3m;
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

  if (lat < 4) {
     c1=4/(3*a*a);
     c2=1/(c*c);
     scon=scoh*(4*pi)*scoh*(4*pi)/(2*a*a*c*sqrt3*econ);
  }
  else if (lat >= 4 and lat <= 5) {
     c1=3/(a*a);
     scon=scoh*(4*pi)*scoh*(4*pi)/(16*a*a*a*econ);
  }
  else if (lat == 6) {
     c1=2/(a*a);
     scon=scoh*(4*pi)*scoh*(4*pi)/(8*a*a*a*econ);
  }
  wint=0;
  t2=hbar/(2*amu*amsc);
  ulim=econ*emax;
  ifl=1;
  nw=maxb;


hexagonalLatticeFactors( a, tsq, c1, c2, lat, tau, nw, tsqx, b, ifl, f, eps, 
  i, wint, twopis, t2, ulim, imax, c );
  // compute lattice factors for hexagonal lattices

  // compute lattice factors for fcc lattices
  //210 continue
  // if (lat > 5) go to 215
  phi=ulim/twopis;
  i1m=int(a*sqrt(phi));
  i1m=15;
  k=0;
   for ( auto i1 = -i1m; i1 < i1m; ++i1 ){
      i2m=i1m;
      for ( auto i2 = -i2m; i2 < i2m; ++i2 ){ 
         i3m=i1m;
         for ( auto i3 = -i3m; i3 < i3m; ++i3 ){
            tsq=taufcc(i1,i2,i3,c1,twothd,twopis);
            if (tsq > 0 and tsq <= ulim) {
               tau=sqrt(tsq);
               w=exp(-tsq*t2*wint)/tau;
               f=w*formf(lat,i1,i2,i3);
               k=k+1;
               if ((2*k) > nw) std::cout << "storage exceeded" << std::endl; 
               b[ifl+2*k-2-1]=tsq;
               b[ifl+2*k-1-1]=f;
            }
          }
        }
   }
   imax=k-1;
   //go to 220


bccLatticeFactors( phi, ulim, twopis, i1m, k, i2m, i3m, tau, w, b, ifl, tsq, 
    f, wint, imax, t2, lat, nw, a, c1 );
   // compute lattice factors for bcc lattices

   // sort lattice factors
  sqrtLatticeFactors( jmin, b, ifl, imax, k, st, sf );
 

  // 220 continue
  for ( auto i = 1; i < imax; ++i ){
      jmin=i+1;
      for ( auto j = jmin; j < k; ++j ){
         if (b[ifl+2*j-2-1] < b[ifl+2*i-2-1]) {
            st=b[ifl+2*i-2-1];
            sf=b[ifl+2*i-1-1];
            b[ifl+2*i-2-1]=b[ifl+2*j-2-1];
            b[ifl+2*i-1-1]=b[ifl+2*j-1-1];
            b[ifl+2*j-2-1]=st;
            b[ifl+2*j-1-1]=sf;
         }
      }
  }
   k=k+1;
   b[ifl+2*k-2-1]=ulim;
   b[ifl+2*k-1-1]=b[ifl+2*k-3-1];
   nw=2*k;

   // convert to practical units and combine duplicate bragg edges.
   bel=-1;
   j=0;
   for ( auto i = 1; i < k; ++i ){
      be=b[ifl+2*i-2-1]*recon;
      bs=b[ifl+2*i-1-1]*scon;
      if (be-bel < toler) {
         b[ifl+2*j-1-1]=b[ifl+2*j-1-1]+bs;
      }
      else {
         j=j+1;
         b[ifl+2*j-2-1]=be;
         b[ifl+2*j-1-1]=bs;
         bel=be;
      }
   }
   nbe=j;
   maxb=2*nbe;
}

