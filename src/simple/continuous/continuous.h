#include <iostream>
#include "simple/generalTools/interpolate.h"
#include "simple/generalTools/integrate.h"
#include "simple/continuous/continuousTools/sumOverTn.h"


template <typename V>
auto getDebyeWaller(const V& P, const V& phononGrid,const int lenRho){
  // Now we want to find lambda_s = int -infty -> infty P(b)*exp(-b/2) db
  // int -inf -> inf P(b)*exp(-b/2)
  // int   0  -> inf P(-b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*2*cosh(b/2) db
  V integrand(lenRho);
  for ( int b = 0; b < lenRho; ++b ){
    integrand[b] = P[b]*2*cosh(phononGrid[b]*0.5);
  }
  return integrate(phononGrid,integrand);
}



template <typename V, typename F>
void continuous(V phononGrid, V rho, V alphas, V betas, const F& kbT, 
  const F& tbeta, const int lat, const int arat){
  int lenRho = rho.size();
  int nAlpha = alphas.size();
  int nBeta  = betas.size();

  // Scale alphas by sc/arat if necessary
  F sc = lat == 1 ? 0.0253/kbT : 1;
  for ( auto& a : alphas ){ a *= sc/arat; }

  // Change phononGrid from eV --> nondimensional beta values    betas = E/kbT
  F inv_kbT = 1.0/(kbT);
  for ( auto& x : phononGrid ){ x *= inv_kbT; }

  // Make sure that rho normalizes to tbeta
  auto invAreaRho = 1.0/integrate(phononGrid, rho);
  for ( auto& y : rho){ y *= tbeta*invAreaRho; }

  V P(phononGrid);
  P[0] = rho[1]/((phononGrid[1]-phononGrid[0])*(phononGrid[1]-phononGrid[0]));
  for ( int b = 1; b < lenRho; ++b ){ 
    P[b] = rho[b]/(2.0*phononGrid[b]*sinh(phononGrid[b]*0.5));
  }

  auto lambda_s = getDebyeWaller(P, phononGrid, lenRho);

  // Calculate T1 = P(b)*exp(-b/2) / lambda_s
  V T1(P.size());
  for ( int b = 0; b < lenRho; ++b ){
    //T1[b] = P[b] * exp(-phononGrid[b]*0.5) / lambda_s;
    T1[b] = P[b] * exp(phononGrid[b]*0.5) / lambda_s;
  }





  /*
  alphas = {2,4};
  betas = {0,1,2,3,4};
  T1 = {1,2,3,2,1};
  lambda_s = 0.5;


  nAlpha = alphas.size();
  nBeta  = betas.size();
  V Tn(nBeta);
  std::copy( T1.begin(), T1.begin() + nBeta, Tn.begin() );

  */
  return;



}



