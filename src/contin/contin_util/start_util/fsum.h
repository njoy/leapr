#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH

template <typename F, typename A> 
double fsum( const int& n, const A& p, const F& tau, const A& betaGrid ){
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
   *   int 0->infty  2 * P(beta) * beta^n * sinh( tau*beta ) dbeta
   *                                or
   *   int 0->infty  2 * P(beta) * beta^n * cosh( tau*beta ) dbeta
   *
   *  depending on whether n is odd of even, respectively.
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * The computed integral is returned
   *
   */

  double func_sum = 0, func_val = 0, dx_left, dx_right;
  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to
                                    // help differentiate betwene sinh and
                                    // cosh while evaluating the integrand
  for( size_t i = 0; i < p.size(); ++i ){
      func_val = even ? 2*p[i]*cosh(betaGrid[i]*tau)*std::pow(betaGrid[i],n) :
                        2*p[i]*sinh(betaGrid[i]*tau)*std::pow(betaGrid[i],n) ;
      dx_left  = 0.0; // We look a half space left, if there's anything <---
      dx_right = 0.0; // We look a half space right, if there's anything --->
      if ( i != 0 )         { dx_left  = (betaGrid[i]-betaGrid[i-1])*0.5; }
      if ( i != p.size()-1 ){ dx_right = (betaGrid[i+1]-betaGrid[i])*0.5; }
      func_val = func_val * ( dx_left + dx_right ); 
    func_sum += func_val;
  } // for i in p
  return func_sum;
}


#endif
