#include <range/v3/all.hpp>
#include "endout/endout_util/inelasticOutput.h"
#include "generalTools/constants.h"
#include "generalTools/tools.h"
#include "ENDFtk.hpp"
#include <iostream>

using namespace njoy::ENDFtk;
using ScatteringFunction   = section::Type<7,4>::TabulatedFunctions::ScatteringFunction;
using ScatteringLaw        = section::Type<7,4>::ScatteringLaw;
using EffectiveTemperature = section::Type<7,4>::EffectiveTemperature;
using Tabulated            = section::Type<7,4>::TabulatedFunctions;
using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;

template <typename Range, typename RangeOfRange >
auto writeInelasticToENDF( const RangeOfRange& fullSAB, const RangeOfRange& fullSAB2,
  const Range alphas, const Range& betas, const Range& temps, const double& za, 
  Range primaryTempf, Range secondaryTempf, int lasym, int lat, int isym, 
  int ilog, ScatteringLawConstants constants, const double& smin ){

  int nbt = (isym == 1 or isym == 3) ? 2*betas.size()-1 : betas.size();
  std::vector<Range> alphaVec (nbt,alphas);
  std::vector< long > boundaries   = { int(alphas.size()) },
                      interpolants = { 4 },  // // ENDF interpolation type 
                      li (temps.size()-1,4);

  std::vector<ScatteringFunction> chunkVectors;
  for (int b = 0; b < nbt; ++b){
    auto boundaries_b   = boundaries;
    auto interpolants_b = interpolants;
    auto li_b           = li;
    auto temps_b        = temps;
    auto alphas_b       = alphas;
    auto out_b          = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,
                                             ilog,lat,b,smin,fullSAB2);
    auto toWrite_b      = std::get<1>(out_b);
    double betaVal = betas[b];

    if ((isym == 1 or isym ==3) and b < betas.size()-1){
        betaVal = -betas[betas.size()-b-1];
    }
    else if ((isym == 1 or isym ==3) and b >= betas.size()-1){
        betaVal = betas[b-betas.size()+1];
    }

    ScatteringFunction chunk_b ( 
                                std::move(betaVal), 
                                std::move(boundaries_b  ),
                                std::move(interpolants_b),
                                std::move(temps_b       ),
                                std::move(li_b          ),
                                std::move(alphas_b      ),
                                std::move(toWrite_b     )
                                );
    chunkVectors.push_back(chunk_b);
  }
  auto awr = constants.atomicWeightRatios()[0];

  ScatteringLaw scatter_law =
    Tabulated( { nbt }, { 4 }, std::move(chunkVectors) );

  auto temps_c = temps;
  auto temps_d = temps;

  EffectiveTemperature principal( { int(temps.size()) }, { 2 }, 
                                  std::move(temps_c),
                                  std::move(primaryTempf)
                                );
  if (constants.numberNonPrincipalScatterers() == 0){
    return section::Type<7,4> ( za, awr, lat, lasym,
                            std::move(constants),
                            std::move(scatter_law),
                            std::move(principal) );
  }
  else {
    if (constants.analyticalFunctionTypes()[0] > 0){
      std::vector< std::optional< EffectiveTemperature > > secondary2 =
                { std::nullopt };

      return section::Type<7,4> ( za, awr, lat, lasym,
                                  std::move(constants),
                                  std::move(scatter_law),
                                  std::move(principal), 
                                  secondary2 );

    }
    else {

      EffectiveTemperature secondary( { int(temps.size()) }, { 2 }, 
                                      std::move(temps_d),
                                      std::move(secondaryTempf)
                                    );
      return section::Type<7,4> ( za, awr, lat, lasym,
                        std::move(constants),
                        std::move(scatter_law),
                        std::move(principal), 
                      { std::optional<EffectiveTemperature>(secondary) } );
    }
  }
}
