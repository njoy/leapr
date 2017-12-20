#include <iostream>
#include <vector>
#include "hexLatticeFactors_util/hexLatticeFactorsHelper.h"
#include "formf.h"
#include "smallFuncs.h"


auto hexLatticeFactors( double a, double tsq, double c1, double c2, 
  int lat, int nw, double tsqx, std::vector<double>& b, int ifl, 
  int i, double wint, double twopis, double t2, double ulim, 
  int imax, double c ){
  double tau, f, eps = 5e-2;
  // compute lattice factors for hexagonal lattices
  int idone;
  double phi=ulim/twopis, w;
  int i1m=int(a*sqrt(phi)), l1, l2, l3, i2m, i3m;
  double w1, w2, w3;
  i1m=i1m+1;
  int k=0;
  for ( auto i1 = 1; i1 <= i1m; ++i1 ){
    l1=i1-1;
    i2m=int((l1+sqrt(3*(a*a*phi-l1*l1)))/2);
    i2m=i2m+1;
    for ( auto i2 = i1; i2 <= i2m; ++i2 ){
      l2=i2-1;
      double x=phi-c1*(l1*l1+l2*l2-l1*l2);
      i3m=0;
      if (x > 0 ) i3m = int(c*sqrt(x));
      i3m=i3m+1;
      for ( auto i3 = 1; i3 <= i3m; ++i3 ){

        l3=i3-1;
        w1=2;
        if (l1 == l2) w1=1;
        w2=2;
        if (l1 == 0 or  l2 == 0) w2=1;
        if (l1 == 0 and l2 == 0) w2=1;
        if (l1 == 0 and l2 == 0) w2=w2/2;
        w3=2;
        if (l3 == 0) w3=1;
        tsq=tausq(l1,l2,l3,c1,c2,twopis);


        if (tsq > 0 and tsq <= ulim) {
          tau=sqrt(tsq);
          w=exp(-tsq*t2*wint)*w1*w2*w3/tau;
          f=w*formf(lat,l1,l2,l3);
          hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, wint, nw, f );
          std::cout << i1 << "    " << i2 << "    " << i3 << "    " << b[13] << "\n" << std::endl;
        }

        tsq = tausq(l1,-l2,l3,c1,c2,twopis);


        if ( i1 == 1 and i2 == 6 and i3 == 2){ return; }


        if (tsq > 0 and tsq <= ulim) {
          tau=sqrt(tsq);
          w=exp(-tsq*t2*wint)*w1*w2*w3/tau;
          f=w*formf(lat,l1,-l2,l3);
          hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, wint, nw, f );
        }
 
       

      } // 3
    } // 2
  } // 1
  imax=k-1;
  //go to 220
}


