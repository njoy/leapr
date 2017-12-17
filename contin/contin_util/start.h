#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>
#include "start_util/start_util.h"

auto start( std::vector<double>& p, double& delta, const double& tev, 
  const double& tbeta ){
  // start makes delta unitless, changing it to what leapr.f90 calls deltab
  
  delta = delta / tev; // make delta is unitless

  // Move phonon distribution rho(beta) to P(beta) by discretely solving at 
  // points delta apart. This follows Eq. 507.
  double beta = delta;
  // What if the first phonon rho value is equal to zero? Then P(beta) is 
  // undefined, which just won't do. So to be safe we approximate the first 
  // value P(beta_0) with a Taylor series. sinh(b/2) ~= b/2 + (b/2)^3/3! + ...,
  // so rho / ( 2 * b * sinh( b/2 ) ) --> rho / ( 2 * b * ( b/2 + ... ) )
  //                                  --> rho / ( b * b )
  
  p[0] = p[1] / ( beta * beta );
  for ( int i = 1; i < p.size(); ++i ){
    p[i] = p[i] / ( 2 * beta * sinh(beta/2) );
    beta += delta; 
  }

  // normalize p so now it integrates to tbeta
  normalize( p, delta, tbeta );

  // calculate debye-waller coefficient and effective temperature
  double lambda_s = fsum( 0, p, 0.5, delta );
  double t_eff    = fsum( 2, p, 0.5, delta ) / ( 2 * tbeta );

  // convert p(beta) --> t1(beta) where t1 is defined to be
  // t1( beta ) = p( beta ) * exp( -beta / 2 ) / lambda_s where
  // lamda_s is the debye-waller coefficient. This relationship
  // is defined by Eq. 525.
  for ( int i = 0; i < p.size(); ++i ){
   p[i] = p[i] * exp( delta * i / 2 ) / lambda_s;
  }
  
  return std::make_tuple( lambda_s, t_eff );
}


