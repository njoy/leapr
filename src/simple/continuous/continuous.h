#include <iostream>
#include "simple/continuous/interpolate.h"
#include "simple/continuous/convolution.h"




template <typename V>
void continuous(V Egrid, V rho, V alphas, V betas){


  // Calculate P(beta) = rho(beta) / beta * sinh(beta/2)
  V P(rho.size());
  P[0] = rho[1] / (Egrid[1]*Egrid[1]);
  for (size_t i = 1; i < rho.size(); ++i ){ 
    P[i] = rho[i]/(Egrid[i]*sinh(Egrid[i]*0.5));
  }


  // Calculate lambda_s (2.0*cosh(b/2) is filling in for e^(-b/2), because 
  // NJOY Eq. 521 is -infty-->infty, but for simplicity we only consider positive
  // Egrid values here (rho is symmetric in energy)
  double lambda_s = 0.0;
  for (size_t i = 0; i < rho.size()-1; ++i ){
    auto val_L = 2.0 * P[i]   * cosh(Egrid[i]*0.5);
    auto val_R = 2.0 * P[i+1] * cosh(Egrid[i+1]*0.5);
    lambda_s += (val_L+val_R)*(Egrid[i+1]-Egrid[i])*0.5;
  }



  // P(beta) to be on the same grid as my requested beta grid 
  V P_beta(rho.size());
  for (size_t b = 0; b < betas.size(); ++b){
    P_beta[b] = interpolate(Egrid,P,betas[b]);
  }


  // calculate T1 (this is now on my beta grid)
  V T1(betas.size());
  std::copy( P_beta.begin(), P_beta.begin() + betas.size(), T1.begin() );
  for (size_t b = 0; b < betas.size(); ++b){
    T1[b] = P_beta[b]*exp(-betas[b]*0.5)/lambda_s;
  }



  V T_left(betas.size());
  std::copy( T1.begin(), T1.begin() + betas.size(), T_left.begin() );
  V T_right = convolve(betas, T1, T_left);

  for (auto& x : T1){ std::cout << x << "   ";} std::cout << std::endl;
  for (auto& x : T_right){ std::cout << x << "   ";} std::cout << std::endl;


  return;
  std::cout << alphas.size() << std::endl;


}




