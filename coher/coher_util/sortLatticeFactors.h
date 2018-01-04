

auto sortLatticeFactors( int ifl, std::vector<double>& b, int k, int jmin, 
  double recon, double& bel, int nbe, int maxb, int j, double& be, double toler, 
  double bs, double scon, double sf, double st, int nw, double ulim, int imax ){

  // Sort lattice factors
  // 220 continue
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
  k = k + 1;
  b[ifl+2*k-2-1] = ulim;
  b[ifl+2*k-1-1] = b[ifl+2*k-3-1];
  nw = 2 * k;
}


