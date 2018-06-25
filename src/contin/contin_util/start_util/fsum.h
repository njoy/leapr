
#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH

#include <range/v3/all.hpp>
#include <iostream>

template<typename T>
T fsum( const int n, const std::vector<T>& p, const T& tau, const T& delta_b ){
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


  T beta = 0, func_sum = 0, func_val = 0;

  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to 
                                    // help differentiate betwene sinh and 
                                    // cosh while evaluating the integrand

  auto b = ranges::view::zip(
             ranges::view::iota(0,int(p.size()))
           | ranges::view::transform([delta_b](auto x){return delta_b*x;}),
	   p );

  /*
  RANGES_FOR(auto entry , b){
    std::cout <<  std::get<0>(entry) << "      " << 
	          std::get<1>(entry) << std::endl;
  }   
  */

  auto funcVal = b | ranges::view::transform([tau,n,even,p,delta_b](auto x){ 
                       T p_i = std::get<1>(x), beta = std::get<0>(x);
		       T value = even ? 
		         2.0 * p_i * cosh(beta*tau) * std::pow(beta,n) :
		         2.0 * p_i * sinh(beta*tau) * std::pow(beta,n) ;
                       return (beta == 0.0 or beta == 1.0*(int(p.size())-1)*delta_b) ?
                         0.5 * value : value;
			 } );



  std::cout << funcVal << std::endl;
  std::cout << std::endl;
  double funcSum = delta_b * ranges::accumulate(funcVal,0.0);
  return funcSum;
  std::cout << "Func Sum: " << funcSum << std::endl;


  if (even) std::cout << "cosh" << std::endl;
  for( size_t i = 0; i < p.size(); ++i ){
    
    func_val = even ? 2.0 * p[i] * cosh( beta * tau ) * std::pow( beta, n ) :
                      2.0 * p[i] * sinh( beta * tau ) * std::pow( beta, n );

    // If at either boundary, cut in half b/c rectangle is only half normal size
    if( i == 0 or i == p.size() -1 ){ func_val = func_val * 0.5 ; }

    std::cout << "---   " << func_val << std::endl; 
    
    beta += delta_b;  func_sum += func_val;

  } // for i in p

  std::cout << func_sum * delta_b << std::endl;
  return func_sum * delta_b;  // return the sum at all requested points, 
                              // multiplied by the width of the rectangles
                              // to give the Riemann summed area
}

#endif
