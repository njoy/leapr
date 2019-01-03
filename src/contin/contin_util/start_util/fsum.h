#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH

template <typename F, typename A> 
double fsum( const int& n, const A& p, const F& tau, const F& delta_b, A betaGrid = A(0) ){
  /* Inputs
   * ------------------------------------------------------------------------
   * n       : appears in equation being evaluated
   * p       : P(beta), which is defined in Eq. 507
   * tau     : appears in equation being evaluated
   * delta_b : increment size used for integral. Doing a Riemann sum with
   *           rectangle sum delta_b wide
   *
   * Operations
   * ------------------------------------------------------------------------
   * * Compute an integral of the form
   *
   *      infty
   *   int       2 * P(beta) * beta^n * sinh( tau*beta ) dbeta
   *      0
   *                                or
   *      infty
   *   int       2 * P(beta) * beta^n * cosh( tau*beta ) dbeta
   *      0
   *
   *  depending on whether n is odd of even, respectively.
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * The computed integral is returned
   *
   */

  bool customGrid = betaGrid.size() > 0 ? true : false;

  double beta = 0, func_sum = 0, func_val = 0, dx_left, dx_right;

  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to
                                    // help differentiate betwene sinh and
                                    // cosh while evaluating the integrand

  for( size_t i = 0; i < p.size(); ++i ){

    if (customGrid){
      func_val = even ? 2 * p[i] * cosh( betaGrid[i]* tau ) * std::pow( betaGrid[i], n ) :
                        2 * p[i] * sinh( betaGrid[i]* tau ) * std::pow( betaGrid[i], n );
      dx_left  = 0.0;
      dx_right = 0.0;
      if ( i != 0 )         { dx_left  = (betaGrid[i]-betaGrid[i-1])/2; }
      if ( i != p.size()-1 ){ dx_right = (betaGrid[i+1]-betaGrid[i])/2; }
      func_val = func_val * ( dx_left + dx_right ); 

    } else {
      func_val = even ? 2 * p[i] * cosh( beta * tau ) * std::pow( beta, n ) :
                        2 * p[i] * sinh( beta * tau ) * std::pow( beta, n );
      beta += delta_b;  
      // If at either boundary, cut in half b/c rectangle is only half normal size
      if( i == 0 or i == p.size() -1 ){ func_val = func_val * 0.5 ; }
    }


    func_sum += func_val;

  } // for i in p

  return customGrid ? func_sum : func_sum * delta_b;
                              // return the sum at all requested points,
                              // multiplied by the width of the rectangles
                              // to give the Riemann summed area
}


#endif
