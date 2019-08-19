#include <iostream> 
#include <vector>
#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>

template<typename Range, typename Float, typename RangeZipped>
auto leapr( int nphon, Float awr, int ncold, Float aws, int lat, Range alpha, Range beta, Range temp_vec, Float delta, Range rho, Float transWeight, Float diffusion_const, Float continWgt, RangeZipped oscEnergiesWeights, Float dka, Range kappaVals ){

  Float lambda_s, tBar, arat=1.0, effectiveTemp, temp, tev, sc, scaling;
  std::vector<Float> sab(alpha.size()*beta.size(),0.0);

  for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
    temp = temp_vec[itemp];
    tev = kb * temp;
    sc = (lat == 1) ? 0.0253/tev : 1.0;
    scaling = sc/arat;

    auto continOutput = 
    contin(nphon, delta, continWgt, scaling, tev, sc, rho, alpha, beta, sab);

    lambda_s = std::get<0>(continOutput);
    tBar     = std::get<1>(continOutput);
    effectiveTemp = temp*tBar;
    

    if (transWeight > 0){
      trans( alpha, beta, transWeight, delta, diffusion_const, sc, scaling, 
        lambda_s, continWgt, effectiveTemp, temp,  sab );
    }

    if (oscEnergiesWeights.size() > 0){
      discre( sc, scaling, lambda_s, transWeight, continWgt, alpha, beta, temp, 
              oscEnergiesWeights, effectiveTemp, sab );
    }


  }

  return sab;

  std::cout << nphon << awr << ncold << aws << lat << alpha[0] << beta[0] << temp_vec[0] << delta << rho[0] << transWeight << diffusion_const << continWgt << std::get<0>(oscEnergiesWeights[0]) << dka << kappaVals[0] << std::endl;

}
