#include <iostream> 
//#include <iomanip>
#include <vector>
//#include <cmath>
//#include <tuple>
//#include <string>
#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include "generalTools/constants.h"
//#include <time.h>

template<typename Range, typename Float>
auto leapr( int nphon, Float awr, int ncold, Float aws, int lat, Range alpha, Range beta, 
  Range temp_vec, Float delta, Range rho, Float trans_weight, Float diffusion_const, Float continWgt, 
  Range oscEnergies, Range oscWeights, Float dka, Range kappaVals ){

  Float lambda_s, tBar, arat=1.0;
  Range symSab(alpha.size()*beta.size());

  for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
    Float temp = temp_vec[itemp];
    Float tev = kb * temp;
    Float sc = (lat == 1) ? 0.0253/tev : 1.0;
    Float scaling = sc/arat;

    auto continOutput = 
    contin(nphon, delta, continWgt, scaling, tev, sc, rho, alpha, beta, symSab);

    lambda_s = std::get<0>(continOutput);
    tBar     = std::get<1>(continOutput);


  }

    //auto lambda_s_t_eff = contin( nphon, delta, tbeta, scaling, tev,
    //  sc, rho, alpha, beta, sym_sab, energyGrid);


  return;
  std::cout << nphon << awr << ncold << aws << lat << alpha[0] << beta[0] << temp_vec[0] << delta << rho[0] << trans_weight << diffusion_const << continWgt << oscEnergies[0] << oscWeights[0] << dka << kappaVals[0] << std::endl;
  // arat is for scaling alpha and beta values later
  //int isecs = 0;
  //bool done = false;

  /*
  Eigen::Tensor<F,3> sym_sab_eigen(alpha.size(),beta.size(),int(temp_vec.size()));
  Eigen::Tensor<F,3> sym_sab_2_eigen(alpha.size(),beta.size(),int(temp_vec.size()));
  sym_sab_eigen.setZero();

  V t_eff_vec ( temp_vec.size(), 0.0 );

  while ( not done ){                      // isecs = 0 --> principal scatterer
    if ( isecs >  0 ){ arat = aws / awr; } // isecs > 0 --> sec. scatterer
      
    for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
      F temp = temp_vec[itemp];
      F tev = bk * temp;
      F sc = 1.0;
      if ( lat == 1 ){ sc = therm/tev; }
      F scaling = sc/arat;

      auto lambda_s_t_eff = contin( nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab_eigen, energyGrid);

      //auto lambda_s_t_eff = contin_NEW( itemp, nphon, tbeta, scaling, tev,
      //  sc, rho, alpha, beta, sym_sab_eigen, energyGrid);


      lambda_s = std::get<0>(lambda_s_t_eff);
      t_eff    = std::get<1>(lambda_s_t_eff);

      if ( trans_weight > 0.0 ){
        trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling,
          itemp, lambda_s, tbeta, t_eff_vec, temp_vec, sym_sab_eigen );
      }

      if ( oscEnergies.size() > 0 ){
        std::vector<std::tuple<F,F>> oscEnergiesWeights(oscEnergies.size());
        for ( size_t i = 0; i < oscEnergies.size(); ++i ){
          oscEnergiesWeights[i] = std::make_tuple(oscEnergies[i],oscWeights[i]);
        }
        discre( itemp, sc, scaling, tev, lambda_s, trans_weight, tbeta, alpha,
          beta, temp_vec, oscEnergiesWeights, t_eff_vec, sym_sab_eigen );
      }

      if ( ncold != 0 ){
	bool free = false;
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, t_eff_vec, 
	  scaling, alpha, beta, dka, kappaVals, lat, free, sym_sab_eigen, 
	  sym_sab_2_eigen );
      }

    }

    done = true;
  }

  return std::make_tuple(lambda_s,t_eff,sym_sab_eigen);
  */

}

