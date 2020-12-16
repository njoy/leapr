#include <range/v3/all.hpp>
#include "endout/endout_util/coherentElasticOutput.h"
#include "ENDFtk.hpp"
#include <iostream>


template <typename Range, typename Float=double>
auto writeCohElasticToENDF( const Range& bragg, const Range& primaryDWF, 
  const Range& secondaryDWF, int numSecondaryScatterers, int secondaryScatterType,
  int numEdges, Range temps, const Float tol = 9e-8 ){ //, const Float& za, const Float& awr){

  using namespace njoy::ENDFtk;
  using CoherentElastic = section::Type<7,2>::CoherentElastic;
  auto braggEnergiesAndXS = processCoherentElastic( bragg, primaryDWF, secondaryDWF, 
    numSecondaryScatterers, secondaryScatterType, numEdges, tol, temps );
  auto energies = std::get<0>(braggEnergiesAndXS);
  auto totalSCR = std::get<1>(braggEnergiesAndXS);
  std::vector<long> boundaries { long(energies.size()) };
  std::vector<long> interpolants { 1 }; // desired ENDF interpolation type
  std::vector<long> temperatureInterpolation (temps.size()-1,2);
  CoherentElastic output( std::move(boundaries), 
                          std::move(interpolants),
                          std::move(temps), 
                          std::move(temperatureInterpolation),
                          std::move(energies), 
                          std::move(totalSCR) ) ;

  return output;
}


