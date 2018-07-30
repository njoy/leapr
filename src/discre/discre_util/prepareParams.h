#include <range/v3/all.hpp>
#include <iostream>

template <typename zipRange>
auto prepareParams( 
  const zipRange& oscEnergiesWeights,
  const double& tev, 
  std::vector<double>& betaVals, double& weight, 
  double& tsave, std::vector<double>& ar, std::vector<double>& dist, 
  std::vector<double>& dbw, const double& bk, 
  std::vector<double>& exb, std::vector<double>& betan, 
  const std::vector<double>& beta, const double& sc ){

  std::vector<double> E(oscEnergiesWeights.size()), W(oscEnergiesWeights.size());
  for ( size_t i = 0; i < oscEnergiesWeights.size(); ++i ){
    E[i] = std::get<0>(oscEnergiesWeights[i]);
    W[i] = std::get<1>(oscEnergiesWeights[i]);
  }
  auto EWRange = ranges::view::zip(E,W);
  auto E_range = EWRange | ranges::view::keys;
  auto W_range = EWRange | ranges::view::values;

  auto sinh_cosh_beta = E_range 
                      | ranges::view::transform([i_tev=1.0/tev](auto E){ 
                           auto b = E * i_tev;
                          return std::make_tuple(b,sinh(0.5*b), cosh(0.5*b));});

  auto cumulativeWeights = ranges::view::iota(1,int(W_range.size())+1)
            | ranges::view::transform([W_range](auto i){ 
                return ranges::accumulate(W_range|ranges::view::slice(0,i),0.0);});

  auto ar_dist_dbw_ranges = 
    ranges::view::concat( 
      ranges::view::zip( W_range, E_range, sinh_cosh_beta ) 
    | ranges::view::transform([](auto t){ 
        auto trig = std::get<2>(t);
        double b = std::get<0>(trig); 
        double w = std::get<0>(t),         
               E = std::get<1>(t), 
          sinh_b = std::get<1>(trig), cosh_b = std::get<2>(trig);
        double coth_b = cosh_b/sinh_b;
        double ar = w/(sinh_b*b), dist = 0.5*w*E*coth_b, dbw = w*coth_b/b;
        return std::make_tuple(b, ar,dist,dbw); }),
      ranges::view::iota(int(E_range.size()),50)
    | ranges::view::transform( [](auto){ return std::make_tuple(0.,0.,0.,0.); }));

  auto ar_ranges = 
    ranges::view::concat( 
      ranges::view::zip( W_range, E_range, sinh_cosh_beta ) 
    | ranges::view::transform([](auto t){ 
        auto trig = std::get<2>(t);
        double b = std::get<0>(trig); 
        double w = std::get<0>(t);
        double sinh_b = std::get<1>(trig);
        double ar = w/(sinh_b*b);
        return ar; }),
      ranges::view::iota(int(E_range.size()),50)
    | ranges::view::transform( [](auto){ return 0.0; }));
  //std::cout << cumulativeWeights << std::endl;
  //std::cout << ar_ranges << std::endl;



  //std::cout << cumulativeWeights << std::endl;
  //std::cout << std::get<0>(ar_dist_dbw_ranges[0]) << std::endl;
  //std::cout << std::get<1>(ar_dist_dbw_ranges[0]) << std::endl;
  //std::cout << std::get<2>(ar_dist_dbw_ranges[0]) << std::endl;
  //std::cout << std::get<3>(ar_dist_dbw_ranges[0]) << std::endl;

  auto betan_range = beta | ranges::view:: transform( [sc](auto b){ return b*sc; } );
  auto exb_range = betan_range | ranges::view::transform( [](auto bn){ 
                                   return exp(-bn); } );

  //std::cout << betan_range << std::endl;
  //std::cout << exb_range << std::endl;


  //auto arRange = ranges::view::transform([](auto t){ return 

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
  //std::cout << std::endl << betaVals[0] << std::endl;
  //std::cout << ar[0] << std::endl;
  //std::cout << dist[0] << std::endl;
  //std::cout << dbw[0] << std::endl;

  for ( size_t b = 0; b < betan.size(); ++b ){
    exb[b] = exp( -beta[b]*sc );
    betan[b] = beta[b]*sc;
  } 
  //std::cout << (exb|ranges::view::all) << std::endl;
  //std::cout << (betan|ranges::view::all) << std::endl;

  return std::make_tuple(ar_dist_dbw_ranges,betan_range,exb_range,ar_ranges);
     
  std::cout << cumulativeWeights << std::endl;

}
