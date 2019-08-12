#include <range/v3/all.hpp>

template <typename Float, typename Range>
auto prepareParams( const std::vector<std::tuple<Float,Float>>& oscEnergiesWeights,
  const Float& tev, Range& betaVals, Float& weight, Float& tsave, Range& ar, 
  Range& dist, Range& dbw, const Float& bk, Range& exb, Range& betan ){
  // Set up oscillator parameters
  using std::exp;
  weight = 0.0;
  tsave = 0.0;
  for ( size_t i = 0; i < oscEnergiesWeights.size(); ++i ){
    Float E = std::get<0>(oscEnergiesWeights[i]);
    Float w = std::get<1>(oscEnergiesWeights[i]);

    betaVals[i] = E / tev;
    weight += w;

    ar[i]   = w / ( sinh(0.5*betaVals[i]) * betaVals[i] );
    dist[i] = 0.5 * w * E / tanh(0.5*betaVals[i]);
    dbw[i]  = w / ( tanh(0.5*betaVals[i]) * betaVals[i] );
    tsave  += dist[i] / bk;
  }

  //betan = ranges::view::iota(0,int(betan.size()))
  //     | ranges::view::transform([sc,beta](auto b){return beta[b]*sc;});

  exb = ranges::view::iota(0,int(betan.size()))
      | ranges::view::transform([betan](auto b){return exp(-betan[b]);});

}
