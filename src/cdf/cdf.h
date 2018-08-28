#include <iostream>
#include "leapr.cpp"


template <typename A>
auto calc_total_ssm(A ssm){

  // INT INT S(a',b') db' da'
  
  double sum = 0;
  for ( int a = 0; a < ssm.dimension(0); ++a ){
    for ( int b = 0; b < ssm.dimension(1); ++b ){
      sum += ssm(a,b,0);
    }
  }
  return sum;
}

template <typename A>
auto calc_eq_14_prime(int b, A ssm){

  // g'(b) = INT  S(a',b) da'
  
  double g_prime = 0;
  for ( int a = 0; a < ssm.dimension(0); ++a ){
    g_prime += ssm(a,b,0);
  }
  return g_prime; 
}


template <typename A>
auto calc_eq_14(int b, A ssm){

  //             INT  S(a',b) da'
  // g(b) =   -----------------------
  //          INT INT S(a',b') db' da'
  
  return calc_eq_14_prime(b,ssm)/calc_total_ssm(ssm);
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


template <typename A>
auto cdf_no_leapr( A ssm ){

  int a_size = ssm.dimension(0);
  int b_size = ssm.dimension(1);

  std::vector<double> eq14(b_size),      eq16(b_size);
  std::vector<double> eq14Prime(b_size), eq16Prime(b_size);
  Eigen::Tensor<double,3> eq15(a_size,b_size,1), eq17(a_size,b_size,1);
  for ( int b = 0; b < b_size; ++b ){
    eq14Prime[b] = calc_eq_14_prime(b,ssm);
    eq16Prime[b] = ( b == 0 ) ? eq14Prime[b] : eq16Prime[b-1] + eq14Prime[b];
    eq14[b] = calc_eq_14(b,ssm);
    eq16[b] = calc_eq_16(b,ssm);
    std::cout << eq16[b] << "      " << eq16Prime[b] << "      " << eq16Prime[b]/calc_total_ssm(ssm) << std::endl;
  }
  for ( int b = 0; b < b_size; ++b ){
    for ( int a = 0; a < a_size; ++a ){
      eq15(a,b,0) = calc_eq_15(a,b,ssm);
      eq17(a,b,0) = calc_eq_17(a,b,ssm);
    }
  }
  return std::make_tuple(eq14,eq15,eq16,eq17);

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
    }
  }
  return std::make_tuple(eq14Values,eq15Values,eq16Values,eq17Values);
}





