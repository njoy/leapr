#include <iostream>
#include <cmath>
#include "trans/trans_util/bessel_K1.h" 

template <typename Float, typename Range>
auto getFreeGas(const Float& trans_weight, const Float& alpha_sc, 
  const int& ndmax, const Float& delta, Range& s_free ){
  using std::exp; using std::pow; using std::sqrt;
  Float beta = 0, wal = trans_weight * alpha_sc;
  int j = 0;

  while (true){
    s_free[j] = exp(-pow(wal-beta,2)/(4*wal))/sqrt(4*M_PI*wal);
    beta += delta; 
    j    += 1;
    if ( (j%2 == 1) and 
         (j+1 >= ndmax or 1e-7*s_free[0] >= s_free[j-1]) ){ return j; }
  }
}


template <typename Float, typename Range>
auto getDiffusion( const Float& trans_weight, const Float& alpha_sc, 
  const int& ndmax, const Float& delta, Range& s_diffusion, 
  const Float& diffusion){

  using std::sqrt; using std::exp;
  int j = 0;
  Float beta = 0.0, wda = trans_weight * diffusion * alpha_sc; 
  while (true){
    Float expTerm = sqrt(beta*beta + 4*wda*wda)*sqrt(diffusion*diffusion+0.25);
    s_diffusion[j] = (2*wda/M_PI)*sqrt(diffusion*diffusion + 0.25) * 
                     bessel_K1_gen(expTerm) / sqrt(beta*beta + 4*wda*wda);
    s_diffusion[j] *= (expTerm <= 1) ? exp(2*wda*diffusion + beta*0.5) :
                                       exp(2*wda*diffusion + beta*0.5 - expTerm); 
    beta += delta; 
    j    += 1;

    if ( (j%2 == 1) and 
         (j+1 >= ndmax or 1e-7*s_diffusion[0] >= s_diffusion[j-1]) ){ return j; }

  } // while
}

