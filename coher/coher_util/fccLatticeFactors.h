#include "formf.h" 
#include "smallFuncs.h" 


auto fccLatticeFactors( const double& twothd, const double& twopis, 
  int lat, std::vector<double>& b, int ifl, double w, int nw,
  double t2, int imax, double c1, double wint, double ulim, double a ){
    

  // compute lattice factors for fcc lattices
  double phi=ulim/twopis;
  double tau;
  int i1m=int(a*sqrt(phi));
  i1m=15;
  int i2m, i3m;
  int k=0;
  for ( auto i1 = -i1m; i1 <= i1m; ++i1 ){
    i2m=i1m;
    for ( auto i2 = -i2m; i2 <= i2m; ++i2 ){ 
      i3m=i1m;
      for ( auto i3 = -i3m; i3 <= i3m; ++i3 ){
        double tsq=taufcc(i1,i2,i3,c1,twothd,twopis);
        if (tsq > 0 and tsq <= ulim) {
          tau=sqrt(tsq);
          w=exp(-tsq*t2*wint)/tau;
          double f=w*formf(lat,i1,i2,i3);
          k=k+1;
          if ((2*k) > nw) std::cout << "storage exceeded" << std::endl; 
          b[ifl+2*k-2-1]=tsq;
          b[ifl+2*k-1-1]=f;
        }
      }
    }
  }
  imax=k-1;
}


