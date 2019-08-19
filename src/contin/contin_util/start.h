#include "contin/contin_util/start_util/normalize.h"
#include "contin/contin_util/start_util/getDebyeWaller.h"
#include "contin/contin_util/start_util/getEffectiveTemp.h"
#include <range/v3/all.hpp>
#include <iostream>


template <typename Range, typename Float>
auto start( const Range& rho, const Float& continWgt, const Range& betaGrid ){
  /* Inputs
   * ------------------------------------------------------------------------
   * p     : excitation frequency spectrum, a function of beta. Originally 
   *         provided by user via Card12
   * delta : desired spacing for the continuous spectrum. This amends beta
   *         values and is originally provided via Card11
   * tev   : temperature in eV. Calculated as T * k_b in leapr.cpp
   * continWgt : normalization for the continuous frequency spectrum. Note that 
   *         the freq. spectrum rho should integrate to 1 (Eq. 508). Well continWgt
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

  // Move phonon distribution rho(beta) to P(beta) by discretely solving at 
  // points delta apart. This follows Eq. 507.
  // What if the first phonon rho value is equal to zero? Then P(beta) is 
  // undefined, which just won't do. So to be safe we approximate the first 
  // value P(beta_0) with a Taylor series. sinh(b/2) ~= b/2 + (b/2)^3/3! + ...,
  // so rho / ( 2 * b * sinh( b/2 ) ) --> rho / ( 2 * b * ( b/2 + ... ) )
  //                                  --> rho / ( b * b )
  // This is is following what the manual instructs on pg. 651 near bottom, that
  // the solid-type spectrum must vary as b^2 as b goes to zero.
  
  auto rhoToP = [](auto betaRho){ 
    auto beta = std::get<0>(betaRho);
    auto rho  = std::get<1>(betaRho);
    return rho / ( 2.0 * beta * sinh(beta*0.5) );
  };
  
  using std::pow;
  auto P = ranges::view::concat(ranges::view::single(rho[1]/pow(betaGrid[1],2)),
                                ranges::view::zip(betaGrid,rho) 
                              | ranges::view::transform(rhoToP) 
                              | ranges::view::drop(1));
  auto betaPZipped_unnormalized = ranges::view::zip(betaGrid, P);
  auto betaPZipped = ranges::view::zip(
                       betaGrid, normalize(betaPZipped_unnormalized, continWgt) );

  Float lambda_s = getDebyeWaller( betaPZipped );
  Float t_eff = getEffectiveTemp(betaPZipped)/(2.0*continWgt);

  auto T1 = betaPZipped 
          | ranges::view::transform( [lambda_s](auto pair){
              auto beta = std::get<0>(pair); 
              auto P = std::get<1>(pair);
              return P * exp(beta*0.5)/lambda_s; });

  return std::make_tuple( lambda_s, t_eff, T1 );
}


