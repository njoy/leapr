#ifndef LEAPR_COHER
#define LEAPR_COHER

#include <iostream>
#include <vector>
#include "generalTools/constants.h"
#include "coher_util/hexLatticeFactors.h"
#include "coher_util/bccLatticeFactors.h"
#include "coher_util/fccLatticeFactors.h"
#include "coher_util/end.h"

auto coher( int iel, int npr, std::vector<double>& b, double maxEnergy ){
  /* Compute Bragg energies and associated structure factors
   * for coherent elastic scattering from graphite, Be, or BeO.
   *
   * inputs :iel tells you the material composition, npr tell you the number of
   * primary atoms, the b vector is the bragg edges vector, nbe is the number 
   * of edges, and maxEnergy (emax) are also there
   */
  int k,imax;
  std::vector<double> 
  //       Graphite   Be         BeO       Al       Pb       Fe  
    aVals {2.4573e-8, 2.2856e-8, 2.695e-8, 4.04e-8, 4.94e-8, 2.8600e-8 },
    cVals {6.7000e-8, 3.5832e-8, 4.39e-8},
    mass  {12.011,    9.0100,    12.50,    26.7495, 207.,    55.454    },
    xsCoh {5.5000,    7.5300,    1.000,    1.49500, 1.00,    12.900    };

  double econ,a=0,c=0,sigmaC,scon,maxTauSq,massScatterer,toler=1.e-6;


  // Lattice Constants (a) in angstroms
  // To get these values, honestly the wikipedia page isn't a bad palce to
  // start. https://en.wikipedia.org/wiki/Lattice_constant. It has all
  // these materials except for Be and BeO. Beryllium can be found at
  // https://www.webelements.com/beryllium/crystal_structure.html, and
  // BeO can be found in the paper "The lattice parameter and density of 
  // beryllium oxide determined by precise x-ray methods" by Bellamy,
  // Baker, and Livel, publlished in Journal of Nuclear Materials in 1962

  // Lattice Constants (c) in angstroms
  // See sources for lattice constants (a) for whereto get these

  // Mass of materials (amu)
  // BeO mass is half off, Al and Fe are a bit off

  // Effective bound coherent scattering cross section (b)
  // To get this, check out "Neutron Scattering Lengths and Cross Sections"
  // by Varley F. Sears, published 1992 in Neutron News.
  // These values are all a little bit off, except for aluminum which is
  // exactly correct, and lead which is comically wrong
  // Lead is way off. Should be like 10 b

  maxEnergy *= ev; // convert this to be in Joules

  econ = ev*8*massNeutron/(1e4*hbar*hbar); // the 1e4 is to make hbar have cm units

  // For hexagonal materials the lattice is described by two constants, a and c,
  // which are defined below. Among other things this is used to define the 
  // reciprocal lattice vector lengths, which is done using tausq. The other
  // reciprocal lattice vector lengths require only one constant.
  // Convert lattice factors from angstroms to cm
  a = aVals[iel-1]; 
  sigmaC = xsCoh[iel-1]/npr;
  massScatterer = amu*mass[iel-1];
  scon = sigmaC*16*M_PI*M_PI;

  //wint=0;                        // This makes me suspicious
  //ifl=1;                         // as does this
  
  maxTauSq = 1e-4*8*massNeutron/(hbar*hbar)*maxEnergy; 
  // max(tau^2) = 8 * m_n * E_max / hbar^2
  //            =     [kg]   [J]  / [J*s]^2
  //            =     [kg] / [kg m^2] = 1/[m^2]
  // The reason there is a 1e-4 here is to make sure that the max tau^2 value
  // is in units of inverse cm^2.

  if ( iel < 4 ){ // compute lattice factors for hexagonal lattices
    c = cVals[iel-1];
    double volume = sqrt(3)*a*a*c/2; // Eq. 559
    scon /= 4*volume*econ;
    imax = hexLatticeFactors( iel, a, c, maxTauSq, b );
  }

  else if ( iel < 6 ){ // compute lattice factors for fcc lattices
    scon/=(16*a*a*a*econ);
    imax = fccLatticeFactors( iel, a, maxTauSq, b );
  } 

  else { // iel == 6  // compute lattice factors for bcc lattices
    scon/=(8*a*a*a*econ);
    imax = bccLatticeFactors( maxTauSq, b, iel, a, massScatterer );
  }
  k = imax + 1;

  int nbe = end( b, k, econ, toler, scon, maxTauSq, imax );
  return std::make_tuple(2*k,2*nbe);
  // first return is the number of nonzero values in b vector
  // second value is 2*#edges

}
#endif
