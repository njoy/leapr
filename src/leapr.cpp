#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include "coher/coher.h"
#include "skold/skold.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>
#include <variant>

template<typename Range, typename Float, typename RangeZipped>
auto leaprTempLoop( int nphon, Float awr, int iel, int npr, int ncold, int lat, Range alpha, Range beta, Range temps, Float delta, Range rho, Float transWgt, Float diffusion_const, Float continWgt, RangeZipped oscEnergiesWgts, Float dka, Range kappaVals, Float cfrac, Float arat ){

  Float tBar, effectiveTemp, temp, tev, sc, scaling;
  std::vector<Float> sab( alpha.size()*beta.size(),0.0),
                     sab2(alpha.size()*beta.size(),0.0);
  Range effectiveTemps(temps.size()), lambdaVals(temps.size());
  std::variant<Range,bool> braggOutput;
  

  for ( size_t itemp = 0; itemp < temps.size(); ++itemp ){ 
    temp = temps[itemp];
    tev = kb * temp;

    sc = (lat == 1) ? 0.0253/tev : 1.0;
    scaling = sc/arat;
    for ( auto& a : alpha ){ a *= scaling; }
    for ( auto& b : beta  ){ b *= sc;      }
    
    delta /= tev;

    //---------------- Incoherent (Elastic and Inelastic) ----------------------
    auto continOutput = contin(nphon, delta, continWgt, rho, alpha, beta, sab);

    lambdaVals[itemp] = std::get<0>(continOutput);
    tBar     = std::get<1>(continOutput);
    effectiveTemp = temp*tBar;
    
    if (transWgt > 0){
      trans( alpha, beta, transWgt, delta, diffusion_const, lambdaVals[itemp], 
             continWgt, effectiveTemp, temp, sab );
    }

    if (oscEnergiesWgts.size() > 0){
      discre( lambdaVals[itemp], transWgt, continWgt, alpha, beta, temp, 
               oscEnergiesWgts, effectiveTemp, sab );
    }

    //--------------- Coherent (Inelastic Approximations) ----------------------
    if (ncold > 0){ 
      bool free = false;
      coldh( tev, ncold, transWgt+continWgt, alpha, beta, dka, kappaVals, free, 
             sab, sab2, tBar );
    }
    else if (kappaVals.size() > 0){
      skold( cfrac, tev, alpha, beta, kappaVals, awr, dka, sab );
    }

    //-------------------------- Coherent (Elastic) ----------------------------
    if (iel > 0){
      double emax = 5.0;
      std::vector<double> bragg ( 60000, 0.0 );
      auto coherOut = coher( iel, npr, bragg, emax );
      int numVals = std::get<0>(coherOut);
      bragg.resize(numVals); 
      braggOutput = bragg;
    }
    else { braggOutput = false; }

    effectiveTemps[itemp] = effectiveTemp;

  }
  return std::make_tuple(sab,effectiveTemps,lambdaVals,braggOutput);
}








template<typename Range, typename Float, typename RangeZipped>
auto leapr( int nphon, Float awr, int iel, int npr, int ncold, Float aws, int lat, Range alpha, Range beta, Range temps, Float delta, Range rho, Float transWgt, Float diffusion_const, Float continWgt, RangeZipped oscEnergiesWgts, Float dka, Range kappaVals, Float cfrac, std::tuple<int,int,Range> secondaryScatterInput ){

  //std::variant<Range,std::tuple<Range,Range>> sab, effectiveTemps, lambdaVals;
  Range sab,  effectiveTemps,  lambdaVals;
  std::variant<Range,bool> braggOutput;

  Range sab2, effectiveTemps2, lambdaVals2;
  std::variant<Range,bool> braggOutput2;

  //-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
  //-/-/-/-/-/-/-/-/-Principal scatterer (isecs = 0)-/-/-/-/-/-/-/-/-/-/-/-/-/-
  //-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
  auto out1 = leaprTempLoop( nphon, awr, iel, npr, ncold, lat, alpha, beta, 
    temps, delta, rho, transWgt, diffusion_const, continWgt, oscEnergiesWgts, 
    dka, kappaVals, cfrac, 1.0 );

  sab             = std::get<0>(out1); effectiveTemps  = std::get<1>(out1);
  lambdaVals      = std::get<2>(out1); braggOutput     = std::get<3>(out1);

  // if number of secondary scatterers is 0 or user wants free gas or diffusion 
  // for the secondary scatterer
  int nss = std::get<0>(secondaryScatterInput),
      b7  = std::get<1>(secondaryScatterInput);
  if (nss == 0 or b7 > 0){
    return std::make_tuple(sab, effectiveTemps, lambdaVals, braggOutput,
                           sab2,effectiveTemps2,lambdaVals2,braggOutput2);
  }

  //-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
  //-/-/-/-/-/-/-/-/-Secondary scatterer (isecs = 1)-/-/-/-/-/-/-/-/-/-/-/-/-/-
  //-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
  
  Range secondaryRho = std::get<2>(secondaryScatterInput);
  auto out2 = leaprTempLoop( nphon, awr, iel, npr, ncold, lat, alpha, beta, 
    temps, delta, secondaryRho, transWgt, diffusion_const, continWgt, 
    oscEnergiesWgts, dka, kappaVals, cfrac, aws/awr );

  sab2            = std::get<0>(out2); effectiveTemps2 = std::get<1>(out2);
  lambdaVals2     = std::get<2>(out2); braggOutput2    = std::get<3>(out2);

  return std::make_tuple(sab, effectiveTemps, lambdaVals, braggOutput,
                         sab2,effectiveTemps2,lambdaVals2,braggOutput2);


}
