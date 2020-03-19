#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"
#include "ENDFtk.hpp"

template <typename Range, typename Float> 
auto processCoherentElastic( const Range& bragg, const Range& dwpix, 
  const Range& dwp1, int numSecondaryScatterers, int secondaryScatterType,
  int numEdges, const Float& tol, const Range& temps ){
  Float w = dwpix[0];
  if (numSecondaryScatterers > 0 and secondaryScatterType == 0){
    w = 0.5*(dwpix[0]+dwp1[0]);
  }

  Range scr (2*numEdges+10,0.0);
  std::vector<Range> totalSCR (temps.size());

  int jmax = 0;
  Float e, sum  = 0.0, suml = 0.0;
  for (int j = 1; j <= numEdges; ++j){
    e = bragg[2*j-2];
    sum += exp(-4.0*w*e)*bragg[2*j-1];
    if ( sum-suml > tol*sum ){
      jmax = j;
      suml = sum;
    }
  }

  Range energies (numEdges+5, 0.0);

  sum = 0;
  for (size_t i = 0; i < temps.size(); ++i){
    sum = 0.0;
    Range scr (numEdges+5,0.0);

    w = dwpix[i];
    if (numSecondaryScatterers > 0 and secondaryScatterType == 0){
      w = 0.5*(dwpix[i]+dwp1[i]);
    }

    int jj = 0;

    if ( i == 0 ){
      for ( int j = 1; j <= numEdges; ++j ){
        if (j <= jmax){ ++jj; }
        e = bragg[2*j-2];
        energies[jj-1] = sigfig(e,7,0);
        sum = sum+exp(-4.0*w*e)*bragg[2*j-1];
        scr[jj-1] = sigfig(sum,7,0);

        if (jj-1 == 56){
            std::cout << bragg[2*j-2] << std::endl;

            //std::cout << energies[56] << "       " << scr[56] << std::endl; 
        }
      }
      energies.resize(jj);
    }
    else {
      for (int j = 1; j <= numEdges; ++j){ 
        if (j <= jmax){ ++jj; }
        e = sigfig(bragg[2*jj-2],7,0);
        sum = sum+exp(-4.0*w*e)*bragg[2*jj-1];
        scr[jj-1] = sigfig(sum,7,0);
      }
    } 
    scr.resize(jj);
    totalSCR[i] = scr;
  }  
  return std::make_tuple(energies,totalSCR);
}


template <typename Range, typename Float>
auto writeCohElasticToENDF( const Range& bragg, const Range& dwpix, 
  const Range& dwp1, int numSecondaryScatterers, int secondaryScatterType,
  int numEdges, const Float& tol, Range temps ){
    //, const Float& za, const Float& awr  ){

  using namespace njoy::ENDFtk;
  using CoherentElastic = section::Type<7,2>::CoherentElastic;
  auto braggEnergiesAndXS = processCoherentElastic( bragg, dwpix, dwp1, 
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

  //section::Type<7,2> 
  //  section(za, awr, CoherentElastic( std::move(boundaries), 
  //                                    std::move(interpolants),
  //                                    std::move(temps), 
  ////                                    std::move(temperatureInterpolation),
  //                                    std::move(energies), 
  //                                    std::move(totalSCR) ) 
  //         );

  //std::string buffer; auto output = std::back_inserter(buffer);
  //section.print(output,27,7); std::cout << buffer << std::endl;


