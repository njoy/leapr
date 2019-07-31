
template <typename Float>
void swap( Float& a, Float& b ){
  Float c = a; a = b; b = c;
}

auto sortLatticeFactors( int ifl, std::vector<double>& b, int& k, 
  int& nw, double ulim, int imax ){

  // Sort lattice factors
  for ( auto i = 1; i <= imax; ++i ){
    for ( auto j = i+1; j <= k; ++j ){
      if (b[ifl+2*j-3] < b[ifl+2*i-3]) {
	swap(b[ifl+2*i-3],b[ifl+2*j-3]);
	swap(b[ifl+2*i-2],b[ifl+2*j-2]);
      }
    }
  }
  k++;
  b[ifl+2*k-3] = ulim;
  b[ifl+2*k-2] = b[ifl+2*k-4];
  nw = 2 * k;
}


auto end( int ifl, std::vector<double>& b, int k, double recon, int& maxb, 
  double toler, double scon, int& nw, double ulim, int imax ){

  std::cout << "NOW IN END" << std::endl;
  std::cout << "  " << std::endl;
  sortLatticeFactors( ifl, b, k, nw, ulim, imax );
  std::cout << "  " << std::endl;
  //std::cout << b[0] << "  " << b[1] << "  " << b[2] << std::endl;
  //std::cout << b[3] << "  " << b[4] << "  " << b[5] << std::endl;
  //std::cout << b[6] << "  " << b[7] << "  " << b[8] << std::endl;
  //std::cout << b[9] << "  " << b[10] << "  " << b[11] << std::endl;
  //std::cout << "  " << std::endl;

  // convert to practical units and combine duplicate bragg edges.

  double bel = -1, be, bs;
  int nbe, j = 0;

  std::cout << recon << "   "  << scon << std::endl;
  for ( auto i = 1; i <= k; ++i ){
    be = b[ifl+2*i-3] * recon;
    bs = b[ifl+2*i-2] * scon;
    if (be-bel < toler) {
      b[ifl+2*j-2] = b[ifl+2*j-2] + bs;
    }
    else {
      j = j+1;
      b[ifl+2*j-3] = be;
      b[ifl+2*j-2] = bs;
      bel=be;
    }
  }
  nbe = j;
  maxb = 2 * nbe;
}


