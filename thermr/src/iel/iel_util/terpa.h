#ifndef THERMR_TERPA_HH
#define THERMR_TERPA_HH
#include "iel/iel_util/terp1.h"
#include <iostream>

auto terpa(double y, double x, std::vector<double> a, int ip, int ir){
 /*--------------------------------------------------------------------
  * Interpolate for y(x) in the TAB1 record packed in a.  Return  0 
  * if x is outside the range of the table.   Here xnext is the next
  * data grid point greater than x.  On entry, ip and ir are starting
  * estimates for the first data point greater than x and the
  * corresponding interpolation range.  Initialize them to 2 and 1
  * before first call to routine.
  *--------------------------------------------------------------------
  */
  int jr, jp, intVar, it, idis;
  double shade=1.00001, xbig = 1.0e12, xnext;
  std::tuple<double,int> out;

  // set up limits and pointers
  ir = std::min(ir, int(round(a[4])) );
  ip = std::min(ip, int(round(a[5])) );
  jr = 5 + 2 * ir;
  jp = 5 + 2 * ( round(a[4]) + ip );
  idis = 0;

  while (true) {
    // locate interpolation interval and law for x
    std::cout << "110  " << ip << "    " << jp << "    " << a[jp-1] << std::endl;
    if (x < a[jp-1]){
      std::cout << "120" << std::endl;
      if (x > a[jp-3]){
        std::cout << "130" << std::endl;
        // interpolate for y in this interval
        intVar=round(a[jr]);
        y = terp1(a[jp-3],a[jp-2],a[jp-1],a[jp],x,intVar);
        if ( intVar == 1 or ( ip != round(a[5]) and a[jp+1] == a[jp-1]) ) {idis=1;}
        out = {a[jp-1],idis};
        return std::tuple<double,int> { a[jp-1], idis };
      }
      if (x == a[jp-3]) {
        std::cout << "140" << std::endl;
        y = a[jp-2];
        intVar = round(a[jr]);
        xnext = a[jp-1];
        if ( intVar == 1 or ( ip != round(a[5]) and a[jp+1] == xnext ) ) {idis=1;}
        return std::tuple<double,int> { xnext, idis };
      }
      if (ip == 2) {
        // special branch for x below first point
        std::cout << "170" << std::endl;
        y = 0;
        return std::tuple<double,int> { a[jp-3], 1 };
      }
      // move down
      jp = jp - 2;
      ip = ip - 1;
      if ( ir != 1 ){
        it=round(a[jr-3]);
        if ( ip <= it ){
          jr = jr-2;
          ir = ir-1;
        }
      }
      //go to 110
    }

    if (ip ==  round(a[5]) ) {
      std::cout << "150" << std::endl;
      // special branch for last point and above
      if (x < shade*a[jp-1]) {
         std::cout << "160" << std::endl;
         y = a[jp];
         xnext = ( y > 0 ) ? shade * shade * a[jp-1] : xbig;
         return std::tuple<double,int> { xnext, idis };
      }
      y = 0;
      return std::tuple<double,int> { xbig, idis };
    }

    // move up
    jp += 2;
    ip += 1;
    it = round(a[jr-1]);
    if ( ip < it ){
      jr += 2;
      ir += 1;
    }
    //go to 110
  }
}

#endif
