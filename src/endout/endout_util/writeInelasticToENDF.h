#include <range/v3/all.hpp>
#include "endout/endout_util/inelasticOutput.h"
#include "generalTools/constants.h"
#include "generalTools/sigfig.h"
#include "ENDFtk.hpp"
#include <iostream>

template <typename Float, typename Range, typename RangeOfRange, 
          typename ScatteringLawConstants >
auto writeInelasticToENDF( const RangeOfRange& fullSAB, const Range alphas, 
  const Range& betas, const Range& temps, const Float& za, 
  Range primaryTempf, Range secondaryTempf,
  int lasym, int lat, int isym, int ilog, ScatteringLawConstants constants ){

  using namespace njoy::ENDFtk;
  using ScatteringFunction   = section::Type<7,4>::Tabulated::ScatteringFunction;
  using ScatteringLaw        = section::Type<7,4>::ScatteringLaw;
  using EffectiveTemperature = section::Type<7,4>::EffectiveTemperature;
  using Tabulated            = section::Type<7,4>::Tabulated;

  std::vector<Range> alphaVec (betas.size(),alphas);
  std::vector< long > boundaries   = { int(alphas.size()) },
                      interpolants = { 4 },  // // ENDF interpolation type 
                      li (temps.size()-1,4);

  std::vector<ScatteringFunction> chunkVectors;
  for (size_t b = 0; b < betas.size(); ++b){
    auto boundaries_b   = boundaries;
    auto interpolants_b = interpolants;
    auto li_b           = li;
    auto temps_b        = temps;
    auto alphas_b       = alphas;
    auto out_b          = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,
                                             ilog,lat,b);
    auto toWrite_b      = std::get<1>(out_b);
    ScatteringFunction chunk_b ( 
                                betas[b], 
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
    Tabulated( { int(betas.size()) }, { 4 }, std::move(chunkVectors) );

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
