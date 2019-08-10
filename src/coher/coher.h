#include <iostream>
#include <vector>
#include "generalTools/constants.h"
#include "coher_util/hexLatticeFactors.h"
#include "coher_util/bccLatticeFactors.h"
#include "coher_util/fccLatticeFactors.h"
#include "coher_util/end.h"

auto coher( int iel, int npr, std::vector<double>& b, 
  double& maxEnergy ){

  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   *
   * inputs :iel tells you the material composition, npr tell you the number of
   * primary atoms, the b vector is the bragg edges vector, nbe is the number 
   * of edges, and maxEnergy (emax) are also there
   */
  int j,k,imax;
  //                        Graphite Be      BeO    Al       Pb    Fe  
  //std::vector<double> aVals {2.4573, 2.2856, 2.695, 4.04,    4.94, 2.86};
  //std::vector<double> cVals {6.700,  3.5832, 4.39};
  //std::vector<double> mass  {12.011, 9.01,   12.5,  26.7495, 207., 55.454};
  //std::vector<double> cohXs {5.50,   7.53,   1.0,   1.495,   1.,   12.9};



  double econ,a=0,c=0,mass,xsCoh,c1,c2,scon,wint,t2,
    maxTauSq,w1,w2,w3,tau,w,f,x,bel,be,bs,

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
  eps = .05;

  maxEnergy *= ev; // convert this to be in Joules
  // to get rid of uninitialized warnings. don't worry i'll change it soon
  imax = 100;



  // initialize.
  econ = 1e-4*ev*8*massNeutron/(hbar*hbar);
  // recon = 1/econ;
  //tsqx = econ/20;

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
  wint=0;                          // This makes me suspicious
  //ifl=1;                         // as does this
  t2=1e4*hbar/(2*amu*mass);
  double massScatterer = amu*mass;
  maxTauSq = 1e-4*8*massNeutron/(hbar*hbar)*maxEnergy; 
  // max(tau^2) = 8 * m_n * E_max / hbar^2
  //            =     [kg]   [J]  / [J*s]^2
  //            =     [kg] / [kg m^2] = 1/[m^2]
  // The reason there is a 1e-4 here is to make sure that the max tau^2 value
  // is in units of inverse cm^2.

  if ( iel < 4 ){ // compute lattice factors for hexagonal lattices
    double volume = sqrt(3)*a*a*c/2; // Eq. 559
    scon /= 4*volume*econ;
    imax = hexLatticeFactors( iel, a, c, maxTauSq, b );
  }

  else if ( iel < 6 ){ // compute lattice factors for fcc lattices
    scon/=(16*a*a*a*econ);
    imax = fccLatticeFactors( iel, a, maxTauSq, massScatterer, b );
  } 

  else { // iel == 6
    // compute lattice factors for bcc lattices
    scon/=(8*a*a*a*econ);
    imax = bccLatticeFactors( maxTauSq, b, iel, a, massScatterer );
  }
  k = imax + 1;

  // nbe is the number of edges
  int nbe = end( b, k, econ, toler, scon, maxTauSq, imax );
  return std::make_tuple(2*k,2*nbe);
  // first return is the number of nonzero values in b vector
  // second value is 2*#edges

}

