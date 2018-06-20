
void swap( double& a, double& b ){
  double c = a;
  a = b; b = c;
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
  k = k + 1;
  b[ifl+2*k-3] = ulim;
  b[ifl+2*k-2] = b[ifl+2*k-4];
  nw = 2 * k;
}


auto end( int ifl, std::vector<double>& b, int k, double recon, int& maxb, 
  double toler, double scon, int& nw, double ulim, int imax ){

  sortLatticeFactors( ifl, b, k, nw, ulim, imax );

  // convert to practical units and combine duplicate bragg edges.

  double bel = -1, be, bs;
  int nbe, j = 0;

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


