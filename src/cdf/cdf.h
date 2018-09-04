#include <iostream>
#include "leapr.cpp"
#include <range/v3/all.hpp>


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
  for ( int aP = 0; aP < ssm.dimension(0); ++aP ){
    denominator += ssm(aP,b,0);
  } 
  if (denominator < 1e-50){ return 0.0; } 
  return numerator/denominator;
}

template <typename A> 
auto calc_eq_15P(int a, int b, A ssm ){
  double numerator   = ssm(a,b,0);
  return numerator;
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

  int a_size = ssm.dimension(0), b_size = ssm.dimension(1);

  std::vector<double> eq14(b_size),  eq16(b_size),
                      eq14P(b_size), eq16P(b_size);
  Eigen::Tensor<double,3> eq15(a_size,b_size,1),  eq17(a_size,b_size,1), 
                          eq15P(a_size,b_size,1), eq17P(a_size,b_size,1);

  for ( int b = 0; b < b_size; ++b ){
    eq14P[b] = calc_eq_14_prime(b,ssm);
    eq16P[b] = ( b == 0 ) ? eq14P[b] : eq14P[b] + eq16P[b-1];
  }

  double inv_T_14_16 = 1.0 / eq16P[b_size-1];

  for ( int b = 0; b < b_size; ++b ){
    eq14[b] = eq14P[b] * inv_T_14_16;
    eq16[b] = eq16P[b] * inv_T_14_16;

    for ( int a = 0; a < a_size; ++a ){
      eq15P(a,b,0) = calc_eq_15P(a,b,ssm);
      eq17P(a,b,0) = ( a == 0 ) ? eq15P(a,b,0) : eq15P(a,b,0) + eq17P(a-1,b,0);
    }
    double inv_T_15_17 = 1.0 / eq17P(a_size-1,b,0);
    for ( int a = 0; a < a_size; ++a ){
      eq15(a,b,0) = eq15P(a,b,0) * inv_T_15_17;
      eq17(a,b,0) = eq17P(a,b,0) * inv_T_15_17;
    }

  }

  std::cout << (eq14|ranges::view::all) << std::endl;
  std::cout << (eq16|ranges::view::all) << std::endl;

  for ( int b = 0; b < b_size; ++b ){
    eq14P[b] = calc_eq_14_prime(b,ssm);
    eq16P[b] = ( b == 0 ) ? eq14P[b] : eq14P[b] + eq16P[b-1];
  }
  auto eq14_rangeP = ranges::view::iota(0,b_size) | 
                     ranges::view::transform([&ssm](auto b){ 
                       return calc_eq_14_prime(b,ssm); } );
  auto eq14_range = eq14_rangeP | 
                    ranges::view::transform([inv_T=1.0/eq16P[b_size-1]]
                      (auto entry){ return entry*inv_T; });
  auto eq16_range = eq16P | 
                    ranges::view::transform([inv_T=1.0/eq16P[b_size-1]]
                      (auto entry){ return entry*inv_T; });

  std::cout << (eq14_range|ranges::view::all) << std::endl;
  std::cout << (eq16_range|ranges::view::all) << std::endl;

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





