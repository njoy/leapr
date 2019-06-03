
template <typename V, typename F>
auto prepareParams( 
  const std::vector<std::tuple<double,double>>& oscEnergiesWeights,
  const F& tev, V& betaVals, F& weight, F& tsave, V& ar, V& dist, V& dbw, 
  const F& bk, V& exb, V& betan, const V& beta, const F& sc ){

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
}
