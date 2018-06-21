
auto prepareParams( 
  const std::vector<std::tuple<double,double>>& oscEnergiesWeights,
 // const std::vector<double>& energy, 
 // const std::vector<double>& weights, 
  const double& tev, 
  std::vector<double>& betaVals, double& weight, 
  double& tsave, std::vector<double>& ar, std::vector<double>& dist, 
  std::vector<double>& dbw, const double& bk, 
  std::vector<double>& exb, std::vector<double>& betan, 
  const std::vector<double>& beta, const double& sc ){

  //std::vector<std::tuple<double,double>> oscEnergiesWeights(energy.size());
  //for ( size_t i = 0; i < energy.size(); ++i ){
  //  oscEnergiesWeights[i] = std::make_tuple(energy[i],weights[i]);
  //}
  
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
