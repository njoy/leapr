#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream> 
#include <iomanip>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include <time.h>


template<typename V, typename F>
auto leapr( int nphon, 
  F awr, int ncold, F aws, 
  int lat, V alpha, V beta, V temp_vec, F delta, V rho, F trans_weight, 
  F diffusion_const, F tbeta, V oscEnergies, V oscWeights, 
  F dka, V kappaVals, V energyGrid = V(0) ){

  F bk    = 8.617385e-5,
    therm = 0.0253,
    lambda_s,t_eff,
    arat  = 1.0; // This is for scaling alpha and beta values later
  // Loop over scatterers and temperatures
  int isecs = 0;
  bool done = false;

  Eigen::Tensor<F,3> sym_sab_eigen(alpha.size(),beta.size(),int(temp_vec.size()));
  Eigen::Tensor<F,3> sym_sab_2_eigen(alpha.size(),beta.size(),int(temp_vec.size()));
  sym_sab_eigen.setZero();

  V t_eff_vec ( temp_vec.size(), 0.0 );

  while ( not done ){
    //if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
    if ( isecs >  0 ){ 
      //std::cout << "Secondary scatterer" << std::endl;
      arat = aws / awr; 
    }
      
    for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
      F temp = temp_vec[itemp];
      F tev = bk * temp;
      F sc = 1.0;
      if ( lat == 1 ){ sc = therm/tev; }
      F scaling = sc/arat;
      if ( itemp == 1 or temp >= 0 ){
       // std::cout << "we want to read in tempdependent parameters" << std::endl;
      } // if 1st temp or some positive temp, we want to calculate
        // the temperature dependent parameters for this specifically 

      auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab_eigen, energyGrid);


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

}

