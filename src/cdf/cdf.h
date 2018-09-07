#include <iostream>
#include "leapr.cpp"
#include <range/v3/all.hpp>


template <typename A>
auto calc14Prime(int b, A ssm){
  // g'(b) = INT  S(a',b) da'
  double g_prime = 0;
  for ( int a = 0; a < ssm.dimension(0); ++a ){ g_prime += ssm(a,b,0); }
  return g_prime; 
}


template <typename A> 
auto calc15Prime(int a, int b, A ssm ){ return ssm(a,b,0); }


template <typename A>
auto cdf_no_leapr( A ssm, double tol=1.0e-20 ){

  int aSize = ssm.dimension(0), bSize = ssm.dimension(1);

  std::vector<double> eq16(bSize);
  A eq17(aSize,bSize,1);

  for ( int b = 0; b < bSize; ++b ){
    eq16[b] = (b == 0) ? calc14Prime(b,ssm) : calc14Prime(b,ssm) + eq16[b-1];
  }
  
//  std::cout << calc14Prime(65,ssm) << "    " << calc14Prime(66,ssm) << "    " 
//            << calc14Prime(67,ssm) << std::endl;
  std::cout << eq16[65] << "   " << eq16[66] << "   "  << eq16[67] << std::endl;

  double inv_T_16 = (eq16[bSize-1] < tol) ? 0.0 : 1.0/eq16[bSize-1];

  for ( int b = 0; b < bSize; ++b ){
    eq16[b] *= inv_T_16;

    for ( int a = 0; a < aSize; ++a ){
      eq17(a,b,0) = (a == 0) ? calc15Prime(a,b,ssm) : 
                               calc15Prime(a,b,ssm) + eq17(a-1,b,0);
    }
    double inv_T_17 = (eq17(aSize-1,b,0) < tol) ? 0.0 : 1.0/eq17(aSize-1,b,0);

    for ( int a = 0; a < aSize; ++a ){ eq17(a,b,0) *= inv_T_17; }
  }

  return std::make_tuple(eq16,eq17);

}



template <typename F, typename I, typename A>
auto cdf(I ntempr, I nphon, I lat, F delta, F twt, F c, F tbeta, A alpha, 
  A beta, A temp, A rho ){

  auto out = leaprWaterSpecific( ntempr, nphon, lat, alpha, beta, temp, delta, 
                                 rho, twt, c, tbeta );

  double lambda_s = std::get<0>(out),
         t_eff    = std::get<1>(out);
  auto eq16 = std::get<3>(out);
  std::cout << eq16[65] << "   " << eq16[66] << "   " << eq16[67] << std::endl;

  auto ssm = std::get<2>(out);
  auto out2 = cdf_no_leapr(ssm);
  return out2;
  return cdf_no_leapr(ssm);
}





