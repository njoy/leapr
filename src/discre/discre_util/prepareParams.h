#include <range/v3/all.hpp>
#include <iostream>

template <typename Float, typename Range,  typename RangeZipped>
auto prepareParams( const RangeZipped& oscEnergiesWeights, const Float& tev, 
  Range& oscBetas, Range& ar, Range& dist, Range& dbw, Range& exb, Range& betan ){

  using std::exp; using std::get;
  Float invTev = 1.0/tev;

  auto oscEnergies = oscEnergiesWeights | ranges::view::keys;
  auto oscWeights  = oscEnergiesWeights | ranges::view::values;

  oscBetas = oscEnergies | ranges::view::transform([invTev](auto E){return E*invTev;});

  auto oscBetasWeights = ranges::view::zip(oscBetas,oscWeights);

  ar  = oscBetasWeights | ranges::view::transform([](auto pair){
              auto beta = get<0>(pair); auto weight = get<1>(pair);
              return weight/(sinh(0.5*beta)*beta); });

  dbw = oscBetasWeights | ranges::view::transform([](auto pair){
              auto beta = get<0>(pair); auto weight = get<1>(pair);
              return weight/(tanh(0.5*beta)*beta); });

  dist = oscEnergiesWeights | ranges::view::transform([invTev](auto pair){
         auto E = get<0>(pair); auto weight = get<1>(pair);
         return 0.5 * weight * E / tanh(0.5*E*invTev); });

  exb = ranges::view::iota(0,int(betan.size()))
      | ranges::view::transform([betan](auto b){return exp(-betan[b]);});

}
