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
auto leaprWaterSpecific( int ntempr, int nphon, int lat, 
  V alpha, V beta, V temp_vec, F delta, V rho, F trans_weight, F diffusion_const, 
  F tbeta ){

  F bk    = 8.617385e-5,
    therm = 0.0253,
    lambda_s,t_eff,
    arat  = 1.0; // This is for scaling alpha and beta values later
  // Loop over scatterers and temperatures
  bool done = false;

  Eigen::Tensor<F,3> sym_sab_eigen(alpha.size(),beta.size(),ntempr);
  Eigen::Tensor<F,3> sym_sab_2_eigen(alpha.size(),beta.size(),ntempr);
  Eigen::Tensor<F,3> eq17(alpha.size(),beta.size(),ntempr);
  sym_sab_eigen.setZero();

  V t_eff_vec ( temp_vec.size(), 0.0 );
  V eq16(beta.size());

  while ( not done ){
      
    for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
      F temp = temp_vec[itemp];
      F tev = bk * temp;
      F sc = 1.0;
      if ( lat == 1 ){ sc = therm/tev; }
      F scaling = sc/arat;

      auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab_eigen);

      lambda_s = std::get<0>(lambda_s_t_eff);
      t_eff    = std::get<1>(lambda_s_t_eff);
      eq16     = std::get<2>(lambda_s_t_eff);

      if ( trans_weight > 0.0 ){
        eq16 = trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling,
          itemp, lambda_s, tbeta, t_eff_vec, temp_vec, sym_sab_eigen );
      }

    }
    done = true;
  }

  return std::make_tuple(lambda_s,t_eff,sym_sab_eigen,eq16,eq17);
  std::cout << diffusion_const << std::endl;
}


template<typename V, typename F>
auto leapr( int nout, std::string title, int ntempr, int iprint, int nphon, 
  int mat, F za, F awr, F spr, int npr, int iel, int ncold, int nss, F aws, 
  int lat, V alpha, V beta, V temp_vec, F delta, int ni, V rho, F trans_weight, 
  F diffusion_const, F tbeta, int nd, V oscEnergies, V oscWeights, int nka, 
  F dka, V kappaVals ){

  F bk    = 8.617385e-5,
    therm = 0.0253,
    lambda_s,t_eff,
    arat  = 1.0; // This is for scaling alpha and beta values later
  // Loop over scatterers and temperatures
  int isecs = 0;
  bool done = false;

  Eigen::Tensor<F,3> sym_sab_eigen(alpha.size(),beta.size(),ntempr);
  Eigen::Tensor<F,3> sym_sab_2_eigen(alpha.size(),beta.size(),ntempr);
  sym_sab_eigen.setZero();

  V t_eff_vec ( temp_vec.size(), 0.0 );

  while ( not done ){
    //if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
    if ( isecs >  0 ){ 
      std::cout << "Secondary scatterer" << std::endl;
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

      std::cout << "going into contin" << std::endl;
      auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab_eigen);

      std::cout << "got back to leapr" << std::endl;

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

  //std::cout << diffusion_const << std::endl;

  std::cout << "Card1: " << nout << std::endl;
  std::cout << "Card2: " <<  title<< std::endl;
  std::cout << "Card3: " << ntempr <<  "     " << iprint << "     " << nphon<< std::endl;
  std::cout << "Card4: " << mat    <<  "     " << za << std::endl;
  std::cout << "Card5: " << awr    <<  "     " << spr    << "     " << npr << "      " << iel << "     " << ncold << std::endl;
  std::cout << "Card6: " << nss    <<  "     " << aws << std::endl;
//  std::cout << "Card7: " << nalpha <<  "     " << nbeta  << "     " << lat<<   std::endl;
  std::cout << "Card8: " << std::endl;
  for( auto entry : alpha   ){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card9: " << std::endl;
  for( auto entry : beta    ){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card10: " << std::endl;
  for( auto entry : temp_vec){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card11: " << delta << "     " << ni <<  std::endl;
  std::cout << "Card12: " << std::endl;
  for( auto entry : rho     ){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card13: " << trans_weight << "     " << diffusion_const << "      " << tbeta << std::endl;
  std::cout << "Card14: " << nd << std::endl;
  std::cout << "Card15: " << std::endl;
  for( auto entry : oscEnergies){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card16: " << std::endl;
  for( auto entry : oscWeights){ std::cout << " -- " << entry << std::endl; }
  std::cout << "Card17: " << nka << "     " << dka << std::endl;
  std::cout << "Card18: " << std::endl;
  for( auto entry : kappaVals){ std::cout << " -- " << entry << std::endl; }

  std::cout << "\n\n\n " << std::endl;

  return std::make_tuple(lambda_s,t_eff,sym_sab_eigen);
}

