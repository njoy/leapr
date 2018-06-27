#include "contin/contin_util/start_util/normalize.h"
#include "contin/contin_util/start_util/fsum.h"
#include <range/v3/all.hpp>

template <typename T, typename V>
auto start( const V& p, T& delta, const T& tev, const T& tbeta ){
  /* Inputs
   * ------------------------------------------------------------------------
   * p     : excitation frequency spectrum, a function of beta. Originally 
   *         provided by user via Card12
   * delta : desired spacing for the continuous spectrum. This amends beta
   *         values and is originally provided via Card11
   * tev   : temperature in eV. Calculated as T * k_b in leapr.cpp
   * tbeta : normalization for the continuous frequency spectrum. Note that 
   *         the freq. spectrum rho should integrate to 1 (Eq. 508). Well tbeta
   *         is the scaling away from 1, if the user wants. Provided in Card13
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Raw freq. spectrum rho(beta) [p] is altered to be P(beta) following 
   *   Eq. 507. It is then normalized and scaled. New P(beta) now exists in p.
   * * Debye-Waller coefficient is calculated, along with effective temperature
   * * P(beta) is turned into T1(beta), following Eq. 525. T1(beta) is in p.
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * T1(beta) exists in the vector p
   * * A tuple containing the Debye-Waller coefficient and effective temperature
   *   is returned.
   *
   */

  delta = delta / tev; // make delta is unitless (leapr.f90 calls this deltab)

  // Move phonon distribution rho(beta) to P(beta) by discretely solving at 
  // points delta apart. This follows Eq. 507.
 
  // What if the first phonon rho value is equal to zero? Then P(beta) is 
  // undefined, which just won't do. So to be safe we approximate the first 
  // value P(beta_0) with a Taylor series. sinh(b/2) ~= b/2 + (b/2)^3/3! + ...,
  // so rho / ( 2 * b * sinh( b/2 ) ) --> rho / ( 2 * b * ( b/2 + ... ) )
  //                                  --> rho / ( b * b )
  // This is is following what the manual instructs on pg. 651 near bottom, that
  // the solid-type spectrum must vary as b^2 as b goes to zero.
  //
  // If it also helps, I guess some intuition you could get is that when
  // we get down to small enough frequencies, everything should turn into the
  // harmonic oscillator that we all knew it was meant to be.
  

  // Make beta grid, [ delta, 2*delta, 3*delta, ... ]
  auto beta = ranges::view::iota(1,int(p.size()))
            | ranges::view::transform([delta](auto x){ return x * delta; });


  // P = rho / (2*beta*sinh(beta/2)). Since we're calculating this using the 
  // rho grid and the beta grid, we need to zip those together so that the 
  // entries are accessible at the same time. Then we add on the 0th term
  // onto the front, which is rho[1]/(beta^2) because of reasons above.
  auto P = ranges::view::concat( 
             ranges::view::single(p[1]/(delta*delta)),
             ranges::view::zip( p | ranges::view::slice(1,int(p.size())), beta )
             | ranges::view::transform( [](auto t){ 
                 T p_i = std::get<0>(t), beta = std::get<1>(t);
                 return p_i / ( 2.0 * beta * sinh(beta/2) ); } ));


  
  // normalize p so now it integrates to tbeta
  auto P1 = normalize( P, delta, tbeta );

  // calculate debye-waller coefficient and effective temperature
  double lambda_sP = fsum( 0, P1, 0.5, delta );
  double t_effP    = fsum( 2, P1, 0.5, delta ) / ( 2 * tbeta );

  // convert p(beta) --> t1(beta) where t1 is defined to be
  // t1( beta ) = p( beta ) * exp( -beta / 2 ) / lambda_s where
  // lamda_s is the debye-waller coefficient. This relationship
  // is defined by Eq. 525.
  auto T1 = ranges::view::zip(P1,ranges::view::iota(0,int(P1.size()))) 
          | ranges::view::transform( [delta,lambda_sP](auto t){ 
              T p_i = std::get<0>(t); int i = std::get<1>(t);
              return p_i * exp(delta*i/2) / lambda_sP; 
            } );
 

  return std::make_tuple( lambda_sP, t_effP, T1 );

}


