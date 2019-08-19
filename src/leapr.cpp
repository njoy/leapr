#include <iostream> 
#include <vector>
#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include "skold/skold.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>

template<typename Range, typename Float, typename RangeZipped>
auto leapr( int nphon, Float awr, int ncold, Float aws, int lat, Range alpha, Range beta, Range temps, Float delta, Range rho, Float transWgt, Float diffusion_const, Float continWgt, RangeZipped oscEnergiesWgts, Float dka, Range kappaVals, Float cfrac ){

  Float lambda_s, tBar, arat=1.0, effectiveTemp, temp, tev, sc, scaling;
  std::vector<Float> sab( alpha.size()*beta.size(),0.0),
                     sab2(alpha.size()*beta.size(),0.0);
  Range effectiveTemps(temps.size());

  for ( size_t itemp = 0; itemp < temps.size(); ++itemp ){ 
    temp = temps[itemp];
    tev = kb * temp;

    sc = (lat == 1) ? 0.0253/tev : 1.0;
    scaling = sc/arat;
    for ( auto& a : alpha ){ a *= scaling; }
    for ( auto& b : beta  ){ b *= sc;      }
    
    delta /= tev;

    //std::cout << "Starting!" << std::endl;
    auto continOutput = 
    contin(nphon, delta, continWgt, rho, alpha, beta, sab);
    //std::cout << "Finished contin!" << std::endl;

    lambda_s = std::get<0>(continOutput);
    tBar     = std::get<1>(continOutput);
    effectiveTemp = temp*tBar;
  //std::cout << (sab|ranges::view::all) << std::endl;
    
    if (transWgt > 0){
      trans( alpha, beta, transWgt, delta, diffusion_const,
        lambda_s, continWgt, effectiveTemp, temp,  sab );
      //std::cout << "Finished trans!" << std::endl;
    }

    if (oscEnergiesWgts.size() > 0){
      discre( lambda_s, transWgt, continWgt, alpha, beta, temp, 
              oscEnergiesWgts, effectiveTemp, sab );
      //std::cout << "Finished discre!" << std::endl;
    }

    if (ncold > 0){ 
      bool free = false;
      coldh( tev, ncold, transWgt, continWgt, scaling, alpha, beta, dka, 
             kappaVals, lat, free, sab, sab2, tBar );
      //std::cout << "Finished coldh!" << std::endl;
    }
    else if (kappaVals.size() > 0){
      skold( cfrac, tev, alpha, beta, kappaVals, awr, dka, scaling, sab );
      //std::cout << "Finished skold!" << std::endl;

    }



    effectiveTemps[itemp] = effectiveTemp;



  }

  return std::make_tuple(sab,effectiveTemps,lambda_s);
  std::cout << (sab|ranges::view::all) << std::endl;

  std::cout << nphon << awr << ncold << aws << lat << alpha[0] << beta[0] << temps[0] << delta << rho[0] << transWgt << diffusion_const << continWgt << std::get<0>(oscEnergiesWgts[0]) << dka << kappaVals[0] << std::endl;

}
