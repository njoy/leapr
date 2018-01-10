//#include "formf.h" 

// BCC
double taubcc( int l1, int l2, int l3, double c1 ){
  return c1 * ( l1*l1 + l2*l2 + l3*l3 + l1*l2 + l2*l3 + l1*l3 ) *
	 4 * M_PI * M_PI;

}

auto bccLatticeFactors( const double& ulim,
  std::vector<double>& b, const int& ifl, const double& wint, 
  const double& t2, const int& lat, const double& a, const double& c1 ){

  // compute lattice factors for bcc lattices
  int im, k = 0;
  double tsq;
  im = int(a*sqrt(ulim/(4*M_PI*M_PI)));
  im = 15;

  for ( auto i1 = -im; i1 <= im; ++i1 ){
    for ( auto i2 = -im; i2 <= im; ++i2 ){
      for ( auto i3 = -im; i3 <= im; ++i3 ){

        tsq = taubcc(i1,i2,i3,c1);

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


