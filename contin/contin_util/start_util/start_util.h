#include <vector>
#include <cmath>


double fsum( const int& n, const std::vector<double>& p, const double& tau, 
             const double& delta_b ){
  /* This function is to compute integrals over the phonon frequency
  * The integral should be from 0 to infinity of some function of the form
  *                 2 * P(beta) * beta^n * hyperbolic
  * integrating with dbeta
  * The hyperbolic function should be either cosh or sinh of ( tau * beta ),
  * depending on whether n is even or odd (respectively).
  * 
  * P(beta) = rho(beta) / ( 2 * beta * sinh(beta/2) ) where rho(beta) is
  * the phonon frequency spectrum.
  */

  double beta = 0, func_sum = 0, func_val = 0;

  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to 
                                    // help differentiate betwene sinh and 
                                    // cosh while evaluating the integrand

  for( int i = 0; i < p.size(); ++i ){
    // Evaluate the integrand. Note that the `even' boolean decides whether
    // to represent using a sinh or cosh. 
    
    func_val = even ? 2 * p[i] * cosh( beta * tau ) * std::pow( beta, n ) :
                      2 * p[i] * sinh( beta * tau ) * std::pow( beta, n );

    if( i == 0 or i == p.size() -1 ){ func_val = func_val * 0.5 ; }

    // increment values for next iteration
    beta += delta_b;  func_sum += func_val;

  } // for i in p

  return func_sum * delta_b;  // return the sum at all requested points, 
                              // multiplied by the width of the rectangles
                              // to give the Riemann summed area
}

std::vector<double> normalize( std::vector<double>& p, const double& delta_b, 
  const double& tbeta ){
  /* normalize first approximates the integral
   * -infty to infty of 2 * P( beta ) * sinh( beta/2 ) 
   * and uses this result to normalize the function P( beta ) to tbeta. 
   * tbeta is user-provided with Card13.
   */

  double sum = fsum( 1, p, 0.5, delta_b ) / tbeta; 
  for ( auto& entry : p ){
    entry = entry / sum;
  } // normalizing loop 

  return p;
}




















