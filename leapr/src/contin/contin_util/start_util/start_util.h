
double fsum( const int& n, const std::vector<double>& p, const double& tau, 
             const double& delta_b ){
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


  double beta = 0, func_sum = 0, func_val = 0;

  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to 
                                    // help differentiate betwene sinh and 
                                    // cosh while evaluating the integrand

  for( size_t i = 0; i < p.size(); ++i ){
    
    func_val = even ? 2 * p[i] * cosh( beta * tau ) * std::pow( beta, n ) :
                      2 * p[i] * sinh( beta * tau ) * std::pow( beta, n );

    // If at either boundary, cut in half b/c rectangle is only half normal size
    if( i == 0 or i == p.size() -1 ){ func_val = func_val * 0.5 ; }

    beta += delta_b;  func_sum += func_val;

  } // for i in p

  return func_sum * delta_b;  // return the sum at all requested points, 
                              // multiplied by the width of the rectangles
                              // to give the Riemann summed area
}

void normalize( std::vector<double>& p, const double& delta_b, 
  const double& tbeta ){
  /* Rearranging Eq. 507 to get a definition for rho(beta), this is the 
   * equation that is being normalized to integrate to tbeta. 
   *
   *             rho(beta) = P(beta) * 2 * beta * sinh( beta / 2 )
   *
   * Inputs
   * ------------------------------------------------------------------------
   * p       : the vector to be normalized. This is P(beta), which is 
   *           defined by Eq. 507. 
   * delta_b : spacing used in the Riemann sum for estimating the integral
   * tbeta   : value to which the above equation should integrate
   *
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Integrate (using fsum) the above equation from 0 to infty.
   * * Scale each P(beta) value so that it is normlized to tbeta
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * P(beta) is amended
   */

  double sum = fsum( 1, p, 0.5, delta_b ) / tbeta; 
  for ( auto& entry : p ){
    entry = entry / sum;
  } 
}




















