

auto do260( std::vector<double> x, std::vector<double> y, double arg, int l, 
  int il ){

  // interpolation section
  double sum = 0, p, pk, in, inp;
  for ( int i = 0; i < il; ++i ){
    p = 1;
    pk = 1;
    in = l + i - 1;
    for ( int ip = 0; ip < il; ++ip ){
      if ( ip == i ){ continue; }
      inp = l + ip - 1;
      p  *= arg   - x[inp];
      pk *= x[in] - x[inp];
    }
    sum += p * y[in] / pk;
  } 
  return sum;
}





auto terp( std::vector<double> x, std::vector<double> y, double arg, 
  int il1 ){

  /*-------------------------------------------------------------------
   * This function does Lagrangian interpolation or
   * extrapolation of ilth order on x and y for arg
   * when the tables are either increasing or decreasing
   *      x          x array, independent variable
   *      y          y array, dependent variable
   *      nl         number of entries in tables of x and y
   *      arg        independent variable value
   *      il         order of interpolation
   *-------------------------------------------------------------------
   */

  int il, l, iadd, il2, ilow, ihi, iusel, iuseh, ibeg, iend, last, m;

  il = il1;
  int nl = x.size();

  // 110
  // If order of interpolation is too large, make it equal to the size of the
  // input vectors
  if ( nl <= il ){
    l = 1; il = nl;                    // Make sure nl is not less than il. 
    return do260( x, y, arg, l, il );  // Continue with interpolation
  } 
   
  // 120
  // If order of interpolation is sufficiently small, continue
  il2 = il / 2;
  // check if tables in increasing or decreasing sequence
  if (x[0] <= x[nl-1]) {
    // increasing sequence
    ilow = il2 + 1;
    ihi = nl - il2 - il%2;
    iusel = 1;
    iuseh = nl - il + 1;
    ibeg = ilow + 1;
    iend = ihi - 1;
    last = iend - il2 + 1;
    iadd = 0;
  }
  else {
    // decreasing sequence
    ilow = nl - il2;
    ihi = il2 + il%2 + 1;
    iusel = nl - il + 1;
    iuseh = 1;
    ibeg = ihi + 1;
    iend = ilow - 1;
    last = 2;
    iadd = 1 - il%2;
  }

  // 170 The following is a weird mess of 170

  // 160
  // If arg is approximately equal to lowest value, return lowest known value
  if ( std::abs( arg - x[ilow-1] ) < arg*1e-10 ) { return y[ilow-1]; }

  // 190
  // If arg is approximately equal to highest value, return highest value
  if ( std::abs( arg - x[ihi-1] ) < arg*1e-10 ) { return y[ihi-1]; }

  // This shows up at the end of 120
  // If arg is lower than x[ilow-1], then set l = iusel and interpolate
  if ( arg < x[ilow-1] ){
    l = iusel;
    return do260( x, y, arg, l, il );
  }

  // This also shows up later
  // If arg is higher than x[ihi-1], then set l = iuseh and interpolate
  if ( arg > x[ihi-1] ){
    l = iuseh;
    return do260( x, y, arg, l, il );
  }

  // 200
  // If arg is a reasonable size, continue
  for ( int n = ibeg; n <= iend; ++n ){

    // If we have a decreasing seqeuence, and the order of interpolation is less
    // than the size of the vector
    m = iusel > 1 ? nl - n + 1 : n;

    // 220 
    
    // 240
    // If arg is approximately equal to a value in my x vector, just return its
    // corresponding value in the y vector
    if ( std::abs(x[m-1] - arg) < 1e-10*arg ){ return y[m-1]; }

    // 250
    // If while iterating through my x vector, I pass over arg, then I should 
    // calculate l and interpolate
    if ( x[m-1] > arg ) { 
      l = m - il2 + iadd;
      return do260( x, y, arg, l, il );
    }
  }
  
  // 230
  // If after iterating through x I haven't found my arg value to interpolate,
  // I assume the latest possible value and interpolate
  l = last;
  return do260( x, y, arg, l, il );
   
}
