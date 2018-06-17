

auto do410( int& i, std::vector<double>& x, const double& xm, const int& nl,
  std::vector<std::vector<double>>& y, const std::vector<double>& yt, const int& j ){
  std::cout << 410 << std::endl;
  x[i]   = x[i-1];
  x[i-1] = xm;
  for ( int il = 0; il < nl; ++il ){
    y[il][i]   =  y[il][i-1];
    y[il][i-1] = yt[il];
  } // enddo

  ++i;
}


