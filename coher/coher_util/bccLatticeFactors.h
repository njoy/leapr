#include "formf.h" 
#include "smallFuncs.h" 


auto bccLatticeFactors( const double& ulim, const double& twopis, 
  std::vector<double>& b, const int& ifl, const double& wint, 
  const double& t2, const int& lat, const double& a, const double& c1 ){

  // compute lattice factors for bcc lattices
  int im, k = 0;
  double tsq;
  im = int(a*sqrt(ulim/twopis));
  im = 15;

  for ( auto i1 = -im; i1 <= im; ++i1 ){
    for ( auto i2 = -im; i2 <= im; ++i2 ){
      for ( auto i3 = -im; i3 <= im; ++i3 ){

        tsq = taubcc(i1,i2,i3,c1,twopis);

        if (tsq > 0 and tsq <= ulim) {
          k += 1;
          if ((2*k) > b.size() ) std::cout << "storage exceeded" << std::endl;
          b[ifl+2*k-2-1] = tsq;
          b[ifl+2*k-1-1] = exp(-tsq*t2*wint)*formf(lat,i1,i2,i3)/sqrt(tsq);
        }

      } // i3
    } // i2 
  } // i1

  return k-1;
}


