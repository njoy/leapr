#include <iostream>
#include <vector>


double formf( int lat, int l1, int l2, int l3 ){
  /* Compute form factors for the specified lattice.
   *       lat=1    graphite
   *       lat=2    Be
   *       lat=3    BeO
   *       lat=4,5  fcc lattice (aluminum, lead)
   *       lat=6    bcc lattice (iron)
   */
   int i;
   double e1,e2,e3, c1=7.54e0, c2=4.24e0, c3=11.31e0, pi=3.14159265358979;
   if (lat == 1) {
      //  graphite
      i=l3/2;
      if ((2*i) != l3) {
        return sin(pi*(l1-l2)/3)*(pi*(l1-l2));
      }
      else {
         return (6+10*cos(2*pi*(l1-l2)/3))/4;
      }
   }
   else if (lat == 2) {
      //  beryllium
      return 1+cos(2*pi*(2*l1+4*l2+3*l3)/6);
    }
   else if (lat == 3) {
      //  beryllium oxide
      return (1+cos(2*pi*(2*l1+4*l2+3*l3)/6))*(c1+c2+c3*cos(3*pi*l3/4));
    }
   else if (lat == 4 or lat == 5) {
      // fcc lattices
      e1=2*pi*l1;
      e2=2*pi*(l1+l2);
      e3=2*pi*(l1+l3);
      return std::pow(1+cos(e1)+cos(e2)+cos(e3),2)+std::pow(sin(e1)+sin(e2)+sin(e3),2);
   }
   else {
      // bcc lattices
      e1=2*pi*(l1+l2+l3);
      return std::pow(1+cos(e1),2)+std::pow(sin(e1),2);
   }
}

