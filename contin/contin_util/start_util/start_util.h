#include <vector>
#include <cmath>


double fsum( const int& n, const std::vector<double>& p, const double& tau, const double& delta_b ){
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

    double exp_increment = exp( delta_b * tau * 0.5 );
    double exp_val = 1, beta = 0, func_sum = 0, func_val = 0;
    int   pos_neg = 1 - 2*(n%2); // +1 if even, -1 if odd. This is to 
                                 // help differentiate betwene sinh and 
                                 // cosh while evaluating the integrand
    for( int i = 0; i < p.size(); ++i ){
        // Evaluate the integrand. Note that the pos_neg decides whether
        // to represent using a sinh or cosh. 
        func_val = ( ( p[i] * exp_val * exp_val )  + 
                     ( p[i] * pos_neg / ( exp_val * exp_val ) ) );

        if( n > 0 ){  func_val = func_val * pow(beta,n);  }

        if( i == 0 or i == p.size() -1 ){ func_val = func_val / 2; }

        beta     = beta + delta_b;
        exp_val  = exp_val * exp_increment;
        func_sum = func_sum + func_val;
    }

    return func_sum * delta_b;  // return the sum at all requested points, 
                                // multiplied by the width of the rectangles
                                // to give the Riemann summed area
    
}

std::vector<double> normalize( std::vector<double>& p, const double& delta_b, const double& tbeta ){
    /* normalize first approximates the integral
     * -infty to infty of 2 * P( beta ) * cosh( beta/2 ) 
     * and uses this result to normalize the function P( beta ) to tbeta. 
     * tbeta is user-provided with Card13.
     */
    double sum = fsum( 1, p, 0.5, delta_b ) / tbeta; 
    for( int i = 0; i < p.size(); i++ ){
        p[i] = p[i] / sum;
    }
    return p;
}




















