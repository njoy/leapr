#include "generalTools/trapezoidIntegral.h"

template <typename Range_Zip>
auto getDebyeWaller(const Range_Zip& beta_P){
  // Now we want to find lambda_s = int -infty -> infty P(b)*exp(-b/2) db
  // int -inf -> inf P(b)*exp(-b/2)
  // int   0  -> inf P(-b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*exp(b/2) + P(b)*exp(-b/2) db
  // int   0  -> inf P(b)*2*cosh(b/2) db
  using std::cosh;
  auto integrand = [](auto xy){ 
    return std::get<1>(xy)*2.0*cosh(std::get<0>(xy)*0.5); };
  return trapezoidIntegral(beta_P,integrand);
}


