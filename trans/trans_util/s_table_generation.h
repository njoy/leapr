#include <iostream>
#include <cmath>
#include <vector>
#include "bessel_K1.h" 


auto free_gas_s_table(const double& trans_weight, const double& alpha_sc, 
  const int& ndmax, const double& delta, std::vector<double>& sd, 
  std::vector<double>& ap ){
  double nsd;
  double beta = 0, wal = trans_weight * alpha_sc;
  int j = 0;

  while ( true ){
    sd[j] = exp( -0.25 * (wal-beta) * (wal-beta)/wal ) / 
            ( 2 * sqrt( M_PI * wal ) );
    ap[j] = beta;
    beta += delta; j += 1;
    // nsd must always be odd for use with Simpson's rule.
    if ( (j%2 == 1) and ( j+1 >= ndmax or 1e-7*sd[0] >= sd[j-1] ) ){
      return j;
    } // if break
  } // while
}


auto diffusion_s_table( const double& diffusion, const double& trans_weight, 
  const double& alpha_sc,  const int& ndmax, const double& delta, 
  std::vector<double>& sd, std::vector<double>& ap ){
    double beta = 0;
    int j = 0;
    double wda = trans_weight * diffusion * alpha_sc; 
    double c4 = 4 * wda * wda;
    double c8 = sqrt(diffusion * diffusion + 0.25) * 2 * wda / M_PI;
    double c3 = 2 * wda * diffusion;
    while ( true ){
    double c6 = sqrt(beta*beta+c4);
    double c7 = c6 * sqrt( diffusion * diffusion + 0.25 );
            
    sd[j] = c7 <= 1 ? c8 * bessel_K1_gen(c7)/c6 * exp( c3 + beta/2 ) :
    c8 * bessel_K1_gen(c7)/c6 * exp( c3 + beta/2 - c7 ) ;
    ap[j] = beta;
    beta += delta;
    j += 1;
    if (j%2 == 1 and ( (j+1 >= ndmax ) or (1e-7*sd[0]>=sd[j-1]) ) ){ 
      return j;
    } // if break
  } // while

}

