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
auto calc_eq_15(int a, int b, A ssm ){
  double numerator   = ssm(a,b,0);
  double denominator = 0;
  for ( int aPrime = 0; aPrime < ssm.dimension(0); ++aPrime ){
    denominator += ssm(aPrime,b,0);
  } 
  if (denominator < 1e-50){ return 0.0; } 
  return numerator/denominator;
}


template <typename A>
auto calc_eq_16(int b, A ssm){
  double sum = 0;
  for ( int i = 0; i <= b; ++i ){
    sum += calc_eq_14(i,ssm);
  }
  return sum;
}

template <typename A> 
auto calc_eq_17( int a, int b, A ssm ){
  double sum = 0;
  for ( int i = 0; i <= a; ++i ){
    sum += calc_eq_15(i,b,ssm);
  }
  return sum;
}

template <typename F, typename I, typename A>
auto cdf(I ntempr, I nphon, I lat, F delta, F twt, F c, F tbeta, A alpha, 
  A beta, A temp, A rho ){

  auto out = leaprWaterSpecific( ntempr, nphon, lat, alpha, beta, temp, delta, 
                                 rho, twt, c, tbeta );

  double lambda_s = std::get<0>(out),
         t_eff    = std::get<1>(out);

  auto ssm = std::get<2>(out);

  A eq14Values(beta.size()), eq16Values(beta.size());
  auto eq15Values = ssm;
  auto eq17Values = ssm;

  for ( int b = 0; b < int(beta.size()); ++b ){ 
    eq14Values[b] = calc_eq_14(b,ssm);
    eq16Values[b] = calc_eq_16(b,ssm);
    for ( int a = 0; a < int(alpha.size()); ++a ){
      eq15Values(a,b,0) = calc_eq_15(b,a,ssm);
      eq17Values(a,b,0) = calc_eq_17(b,a,ssm);
      //if ( b == 0 ){ std::cout << eq17Values(a,b,0) << std::endl; }
    }
  }
  return std::make_tuple(eq14Values,eq15Values,eq16Values,eq17Values);
}
