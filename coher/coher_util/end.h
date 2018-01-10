
void swap( double& a, double& b ){
  double c = a;
  a = b; b = c;
}


auto sortLatticeFactors( std::vector<double>& b, int ifl, int imax,
  int k, double& ulim ){
  // sort lattice factors
  int jmin; 
  double st, sf;
  for ( auto i = 1; i < imax; ++i ){
    jmin=i+1;
    for ( auto j = jmin; j < k; ++j ){
      if (b[ifl+2*j-2-1] < b[ifl+2*i-2-1]) {
        st = b[ifl+2*i-2-1];
        sf = b[ifl+2*i-1-1];
        b[ifl+2*i-2-1] = b[ifl+2*j-2-1];
        b[ifl+2*i-1-1] = b[ifl+2*j-1-1];
        b[ifl+2*j-2-1] = st;
        b[ifl+2*j-1-1] = sf;
      }
    }
  }
  k += 1;
  b[ifl+2*k-3] = ulim;
  b[ifl+2*k-1] = b[ifl+2*k-4];
  return 2*k;
}



auto end( int ifl, std::vector<double>& b, int k, double recon, int maxb, 
  double toler, double scon, int nw, double ulim, int imax ){

  sortLatticeFactors( b, ifl, imax, k, ulim );

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


