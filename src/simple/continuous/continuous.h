#include <iostream>
#include "simple/continuous/interpolate.h"
#include "simple/continuous/integrate.h"




template <typename V, typename F>
void continuous(V phononGrid, V rho, V alphas, V betas, F kbT){
  int nAlpha = alphas.size();
  int nBeta  = betas.size();

  // Change phononGrid from eV --> nondimensional beta values    betas = E/kbT
  F inv_kbT = 1.0/(kbT);
  for ( auto& x : phononGrid ){ x *= inv_kbT; }

  // Change rho, which is now on the phononGrid, to be on the betas grid
  V P(nBeta);
  P[0] = rho[1]/(betas[1]*betas[1]);
  for ( int b = 1; b < nBeta; ++b ){ 
    F beta = betas[b];
    P[b] = interpolate(phononGrid,rho,beta)/(2.0*beta*sinh(beta*0.5));
  }

  // Now we want to find lambda_s = int -infty -> infty P(b)*exp(-b/2) db
  // int -inf -> inf P(b)*exp(-b/2)
  // int   0  -> inf P(-b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*2*cosh(b/2) db
  V integrand(P.size());
  for ( size_t b = 0; b < P.size(); ++b ){
    integrand[b] = P[b]*2*cosh(betas[b]*0.5);
  }
  F lambda_s = integrate(betas,integrand);

  // Calculate T1 = P(b)*exp(-b/2) / lambda_s
  V T1(P.size());
  for ( size_t b = 0; b < P.size(); ++b ){
    T1[b] = P[b] * exp(-betas[b]*0.5) / lambda_s;
  }

  //print(betas);
  //print(P);
  //print(T1);


  alphas = {2,4};
  betas = {0,1,2,3,4};
  T1 = {1,2,3,2,1};
  lambda_s = 0.5;


  nAlpha = alphas.size();
  nBeta  = betas.size();
  V Tn(nBeta);
  std::copy( T1.begin(), T1.begin() + nBeta, Tn.begin() );






  return;

  std::cout << alphas.size() << std::endl;


}




