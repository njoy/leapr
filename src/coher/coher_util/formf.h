
#ifndef COHER_COHERUTIL_FORMF
#define COHER_COHERUTIL_FORMF

#include <iostream>
#include <vector>


inline double formf( int lat, int l1, int l2, int l3 ){
  /* Overview
   * ------------------------------------------------------------------------
   * Computes form factors for the specified lattice.
   *       lat = 1    graphite
   *       lat = 2    Be
   *       lat = 3    BeO
   *       lat = 4,5  fcc lattice (aluminum, lead)
   *       lat = 6    bcc lattice (iron)
   *
   * Inputs
   * ------------------------------------------------------------------------
   * lat : Integer indicating the material
   * l1  : Iterator for the x axis
   * l2  : Iterator for the y axis
   * l3  : Iterator for the z axis
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Each of the above structures has a certain number of atoms per unit 
   *   cell, at specific relative positions. The equations below take those
   *   positions, creates corresponding phi values (See Eq. 560-563) and then
   *   solves Eq. 554. Note that these intermediate phi values that are 
   *   implicitly computed are defined to be 
   *                   phi_j = \vec{tau} \cdot \vec{rho_j}
   *   where tau is the vector of one particular "shell" of the reciprocal
   *   lattice, and rho_j are the position vectors of atoms. Thus, phi_j
   *   represetns the phases for the atoms. (See pg. 659)
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * Form factor |F|^2, as defined in Eq. 554.
   */

  using std::pow;

   int i;
   double e1,e2,e3, c1=7.54e0, c2=4.24e0, c3=11.31e0, pi=3.14159265358979;

   //  graphite -- HEX
   if (lat == 1) {
     // This follows Eq. 564. F is set to either value depending on if l3 is 
     // even or odd. Note that the returned value differs from Eq. 564 however
     // in that either possibility is divided by 4.
     return l3%2 == 0 ? (6+10*cos(2*pi*(l1-l2)/3))/4 :  // even
	                pow(sin(pi*(l1-l2)/3),2);  // odd
   }

   //  beryllium -- HCP
   else if (lat == 2) {
     // This follows Eq. 565, since this is a hexagonal close packed (hcp) 
     // structure. Notice that the returned value differs from Eq. 565 in that 
     // it is divided by 2.
     return 1+cos(2*pi*(2*l1+4*l2+3*l3)/6);
   }

   //  beryllium oxide -- 2 Interpenetrating HCP
   else if (lat == 3) {
     // This follows Eq. 566 in spirit, with the slight caveat that the 
     // c1 = r1^2, c2 = r2^2, so c1 + c2 = r^2 and c3 = 2*r1*r2
     return (1+cos(2*pi*(2*l1+4*l2+3*l3)/6))*(c1+c2+c3*cos(3*pi*l3/4));
   }

   // fcc lattices
   else if (lat == 4 or lat == 5) {
     // This doesn't follow a particular equation, as it is left as an exercise
     // to the reader. But it seems as if it's stating that there are 4 atoms
     // per unit cell at (0,0,0) (1,0,0) (1,1,0) (1,0,1)
     e1 = 2 * pi * l1;
     e2 = 2 * pi * (l1+l2);
     e3 = 2 * pi * (l1+l3);
     return pow( cos(e1) + cos(e2) + cos(e3) + 1, 2 ) + 
            pow( sin(e1) + sin(e2) + sin(e3),     2 );
   }

   // bcc lattices
   else {
      return pow( cos( 2*pi*(l1+l2+l3) ) + 1,2 ) +
	     pow( sin( 2*pi*(l1+l2+l3) ),    2 );
   }
}

#endif
