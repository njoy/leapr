#include <iostream>
#include "leapr.cpp"


template <typename A>
auto calc_total_ssm(A ssm){
  double sum = 0;
  for ( int a = 0; a < ssm.dimension(0); ++a ){
    for ( int b = 0; b < ssm.dimension(1); ++b ){
      sum += ssm(a,b,0);
    }
  }
  return sum;
}

template <typename A>
auto calc_eq_14(int b, A ssm){
  double numerator = 0;
  for ( int a = 0; a < ssm.dimension(0); ++a ){
    numerator += ssm(a,b,0);
  }
  return numerator/calc_total_ssm(ssm);
}

template <typename A> 
auto calc_eq_18(int a, int b, A ssm ){
  double numerator   = ssm(a,b,0);
  double denominator = 0;
  for ( int aPrime = 0; aPrime < ssm.dimension(0); ++aPrime ){
    denominator += ssm(aPrime,b,0);
  } 
  return numerator/denominator;
}


template <typename A>
auto calc_eq_16(int b, A ssm){
  double sum = 0;
  for ( int i = 0; i < b; ++i ){
    sum += calc_eq_14(i,ssm);
  }
  return sum;
}

template <typename F, typename I, typename A>
I cdf(I ntempr, I nphon, 
  I lat, 
  F delta, F twt, F c, F tbeta, A alpha, A beta, A temp, A rho
  ){

    auto out = leaprWaterSpecific( ntempr, nphon, lat, alpha, 
        beta, temp, delta, rho, twt, c, tbeta );
 

    double lambda_s = std::get<0>(out),
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    for ( int i = 0; i < int(beta.size()); ++i ){ 
      //std::cout << beta[i] << "           " << calc_eq_14(i,ssm) << std::endl;
      //std::cout << calc_eq_14(i,ssm) << std::endl;
    }
    std::cout << std::setprecision(16) << calc_eq_14(25,ssm) << std::endl;
    std::cout << std::endl;

  return 0;
}
