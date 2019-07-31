#include <iostream>
#include <vector>
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
  int i,j,k,imax,ifl,nw;
  double amne,econ,tsqx,a=0,c=0,mass,xsCoh,c1,c2,recon,scon,wint,t2,
    ulim,phi,w1,w2,w3,tsq,tau,w,f,x,bel,be,bs,

    // Lattice Constants (a) in angstroms
    // To get these values, honestly the wikipedia page isn't a bad palce to
    // start. https://en.wikipedia.org/wiki/Lattice_constant. It has all
    // these materials except for Be and BeO. Beryllium can be found at
    // https://www.webelements.com/beryllium/crystal_structure.html, and
    // BeO can be found in the paper "The lattice parameter and density of 
    // beryllium oxide determined by precise x-ray methods" by Bellamy,
    // Baker, and Livel, publlished in Journal of Nuclear Materials in 1962
    aGraphite = 2.4573,
    aBe       = 2.2856,
    aBeO      = 2.695,
    aAl       = 4.04,
    aPb       = 4.94,
    aFe       = 2.86,

    // Lattice Constants (c) in angstroms
    // See sources for lattice constants (a) for whereto get these
    cGraphite = 6.700,
    cBe       = 3.5832,
    cBeO      = 4.39,

    // Mass of materials (amu)
    massGraphite = 12.011,  
    massBe       = 9.01,   
    massAvgBeO   = 12.5,     // 1/2 BeO
    alMass       = 26.7495,  // (a bit off?)
    massPb       = 207.,        
    massFe       = 55.454,   // (a bit off?)

    // Effective bound coherent scattering cross section (b)
    // To get this, check out "Neutron Scattering Lengths and Cross Sections"
    // by Varley F. Sears, published 1992 in Neutron News.
    // These values are all a little bit off, except for aluminum which is
    // exactly correct, and lead which is comically wrong
    xsCohGraphite = 5.50, 
    xsCohBe       = 7.53,
    xsCohBeO      = 1.0,
    xsCohAl       = 1.495,
    xsCohPb       = 1.,     // Right out. Should be ~10 b
    xsCohFe       = 12.9,   


    toler = 1.e-6,
    eps = .05,
    amassn = 1.008664904,   // mass of neutron in amu
    amu = 1.6605402e-27,
    ev = 1.60217733e-19,
    hbar = 1.05457266e-34;  // [J*s]

  // to get rid of uninitialized warnings. don't worry i'll change it soon
  phi = 0; tsq = 0; i = 0; imax = 100;



  // initialize.
  amne = amassn * amu;  // Mass of neutron in kg ( 1.6749286E-27 )
  econ = 1e-4*ev*8*amne/(hbar*hbar);
  recon = 1/econ;
  tsqx = econ/20;

  // For hexagonal materials the lattice is described by two constants, a and c,
  // which are defined below. Among other things this is used to define the 
  // reciprocal lattice vector lengths, which is done using tausq. The other
  // reciprocal lattice vector lengths require only one constant.
  if (iel == 1){
     // graphite constants.
     a = aGraphite;
     c = cGraphite;
     mass = massGraphite;
     xsCoh = xsCohGraphite / npr;
  }
  else if (iel == 2) {
     //  beryllium constants
     a = aBe;
     c = cBe;
     mass = massBe;
     xsCoh = xsCohBe / npr;
  }
  else if (iel == 3) {
     //  beryllium oxide constants
     a = aBeO;
     c = cBeO;
     mass = massAvgBeO;
     xsCoh = xsCohBeO / npr;
  }
  else if (iel == 4) {
     // aluminum constants
     a = aAl;
     mass = alMass;
     xsCoh = xsCohAl / npr;
  }
  else if (iel == 5) {
     // lead constants
     a = aPb;
     mass = massPb;
     xsCoh = xsCohPb / npr;
  }
  else if (iel == 6) {
     // iron constants
     a = aFe;
     mass = massFe;
     xsCoh = xsCohFe / npr;
  }
  else { 
    std::cout << "oh no, not a valid iel value" << std::endl;
    throw std::exception(); 
  }

  // Convert lattice factors from angstroms to cm
  a *= 1e-8; c *= 1e-8;

  scon = xsCoh*16*M_PI*M_PI;
  if (iel < 4) {
     c1=4/(3*a*a);
     c2=1/(c*c);
     double volume = sqrt(3)*a*a*c/2;
     scon /= 4*volume*econ;
     //scon /= (2*a*a*c*sqrt(3)*econ);
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
  t2=1e4*hbar/(2*amu*mass);
  ulim=econ*emax;
  ifl=1;
  nw=maxb;


  if ( iel < 4 ){
    // compute lattice factors for hexagonal lattices
    phi = ulim/(4.0*M_PI*M_PI);
    int i1m = a*sqrt(phi) + 1;
    imax = hexLatticeFactors( a, tsq, c1, c2, iel, nw, tsqx, b, ifl,  
    i, wint, t2, ulim, c, i1m );
    k = imax + 1;
  }

  else if ( iel < 6 ){
    // compute lattice factors for fcc lattices
    imax = fccLatticeFactors( iel, b, ifl, nw, t2, c1, wint, ulim, a );
    k = imax + 1;
  } 

  else { // iel == 6
    // compute lattice factors for bcc lattices
    imax = bccLatticeFactors( ulim, b, ifl, wint, t2, iel, a, c1 );
    k = imax + 1;
  }

  // nbe is the number of edges
  int nbe = end( ifl, b, k, recon, maxb, toler, scon, nw, ulim, imax );
  maxb = 2*nbe;
 
  //nbe = 1;
  return nbe;

}

