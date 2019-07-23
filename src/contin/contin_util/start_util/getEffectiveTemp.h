#include "generalTools/trapezoidIntegral.h"

template <typename Range_Zip>
auto getEffectiveTemp(const Range_Zip& beta_P){
  using std::cosh; using std::pow;
  auto integrand = [](auto xy){ 
    return std::get<1>(xy)*pow(std::get<0>(xy),2)*2.0*cosh(std::get<0>(xy)*0.5); };
  return trapezoidIntegral(beta_P,integrand);
}


