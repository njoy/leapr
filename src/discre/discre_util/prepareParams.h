#include <iostream>
#include <range/v3/all.hpp>

// Float, Array, Array of Tuples
template <typename F, typename A, typename A_of_T>
auto prepareParams( const A_of_T& oscEnergiesWeights, const F& tev, A& betaVals, 
  F& weight, F& tsave, A& ar, A& dist, A& dbw, const F& bk, A& exb, A& betan, 
  const A& beta, const F& sc ){

  std::vector<double> eVector(oscEnergiesWeights.size());
  std::vector<double> wVector(oscEnergiesWeights.size());
  int counter = 0;
  for ( auto entry : oscEnergiesWeights ){ 
    eVector[counter]=std::get<0>(entry);
    wVector[counter++]=std::get<1>(entry);
  }
  auto E_range = ranges::view::zip(eVector,wVector) | ranges::view::keys;
  auto W_range = ranges::view::zip(eVector,wVector) | ranges::view::values;

  auto betaVals_range = 
    ranges::view::concat(
      E_range | ranges::view::transform([i_tev=1.0/tev](auto E){ 
        return E*i_tev; }),
      ranges::view::iota(int(E_range.size()),50) 
    | ranges::view::transform( [](auto){ return 0.0; } )
    );

  auto sinh_cosh_beta = betaVals_range | ranges::view::transform([](auto b){ 
                          return std::make_tuple(b,sinh(0.5*b), cosh(0.5*b));});

  auto cumulativeWeights = ranges::view::iota(1,int(W_range.size())+1)
            | ranges::view::transform([W_range](auto i){ 
                return ranges::accumulate(W_range|ranges::view::slice(0,i),0.0);});

  auto ar_range = ranges::view::concat( ranges::view::zip( W_range, E_range, sinh_cosh_beta ) 
                | ranges::view::transform([](auto t){ 
                    auto w = std::get<0>(t); 
                    auto b = std::get<0>(std::get<2>(t)); 
                    auto sinh_b = std::get<1>(std::get<2>(t));
                    return w / ( sinh_b * b ); } ), 
                          ranges::view::iota(int(E_range.size()),50)
                        | ranges::view::transform( [](auto){ return 0.0; } ));
  auto dist_range = ranges::view::concat( ranges::view::zip( W_range, E_range, sinh_cosh_beta ) 
                | ranges::view::transform([](auto t){ 
                    auto w = std::get<0>(t); 
                    auto E = std::get<1>(t); 
                    auto sinh_b = std::get<1>(std::get<2>(t));
                    auto cosh_b = std::get<2>(std::get<2>(t));
                    return 0.5 * w * E * cosh_b / sinh_b; } ), 
                          ranges::view::iota(int(E_range.size()),50)
                        | ranges::view::transform( [](auto){ return 0.0; } ));
  auto dbw_range = ranges::view::concat( ranges::view::zip( W_range, E_range, sinh_cosh_beta ) 
                | ranges::view::transform([](auto t){ 
                    auto w = std::get<0>(t); 
                    auto b = std::get<0>(std::get<2>(t)); 
                    auto sinh_b = std::get<1>(std::get<2>(t));
                    auto cosh_b = std::get<2>(std::get<2>(t));
                    return w * cosh_b / ( sinh_b * b ); } ), 
                          ranges::view::iota(int(E_range.size()),50)
                        | ranges::view::transform( [](auto){ return 0.0; } ));

  auto betan_range = beta | ranges::view:: transform( [sc](auto b){ return b*sc; } );
  auto exb_range = betan_range | ranges::view::transform( [](auto bn){ 
                                   return exp(-bn); } );



  //std::cout << exb_range << std::endl;
  // Set up oscillator parameters
  weight = 0.0;
  tsave = 0.0;
  for ( size_t i = 0; i < oscEnergiesWeights.size(); ++i ){
    double E = std::get<0>(oscEnergiesWeights[i]);
    double w = std::get<1>(oscEnergiesWeights[i]);

    betaVals[i] = E / tev;
    weight += w;

    ar[i]   = w / ( sinh(0.5*betaVals[i]) * betaVals[i] );
    dist[i] = 0.5 * w * E / tanh(0.5*betaVals[i]);
    dbw[i]  = w / ( tanh(0.5*betaVals[i]) * betaVals[i] );
    tsave  += dist[i] / bk;
  }

  for ( size_t b = 0; b < betan.size(); ++b ){
    exb[b] = exp( -beta[b]*sc );
    betan[b] = beta[b]*sc;
  } 
  //std::cout << (ar|ranges::view::all) << std::endl;
  //std::cout << (dist|ranges::view::all) << std::endl;
  //std::cout << (dbw|ranges::view::all) << std::endl;
}
