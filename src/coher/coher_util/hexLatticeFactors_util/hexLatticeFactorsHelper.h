
template <typename Float, typename Range>
auto hexLatticeFactorsHelper( int& k, const Float& tsq, const Float& tsqx, 
  Range& b, const int& ifl, const int& nw, const Float& f, int& i ){
 
  if (k <= 0 or tsq <= tsqx) {
    k += 1;
    if ((2*k) > nw) { throw std::invalid_argument("2k must be <= nw"); } 
    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;
    return;
  }

  else {
    for ( i = 1; i < k+1; ++i ){
      if ( tsq >= b[ifl+2*i-3] and tsq <= 1.05 * b[ifl+2*i-3] ) {
        b[ifl+2*i-2] += f;
	return;
      } // if
    } // while
    k++;
    b[ifl+2*k-3] = tsq;
    b[ifl+2*k-2] = f;
  } // else
}
