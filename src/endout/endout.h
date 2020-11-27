#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"
#include "coherentElastic/coherentElastic.h"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "endout/endout_util/writeIncElasticToENDF.h"
#include "endout/endout_util/writeCohElasticToENDF.h"
#include "ENDFtk.hpp"
#include "lipservice.hpp"

using namespace njoy::ENDFtk;
using Elastic         = section::Type<7,2>;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using Inelastic       = section::Type<7,4>;
using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;

template <typename Range> 
auto scaleDebyeWallerCoefficients( int numSecondaryScatterers, 
  int secondaryScatterType, Range& dwpix, Range& dwp1, const Range& temps, 
  const Range& awrVec ){
  // display endf t-effective and debye-waller integral
  auto awr = awrVec[0],
       aws = awrVec[1];
  for (size_t i = 0; i < temps.size(); ++i){
    if (numSecondaryScatterers == 0 or secondaryScatterType > 0){
       dwpix[i] /= (awr*temps[i]*kb);
    }
    else {
       dwpix[i] /= (aws*temps[i]*kb);
       dwp1[i]  /= (awr*temps[i]*kb);
    }
  }
}



template <typename Range>
auto endoutNew( const nlohmann::json& jsonInput,
  std::vector<Range>& sab,
  const std::vector<Range>& principalScatterSAB, const Range& alphas, 
  const Range& betas, Range& dwpix, Range& dwp1, 
  const Range& bragg, int numEdges, 
  Range primaryTempf, Range secondaryTempf ){ 
  using std::pow;


  unsigned int npr = jsonInput["npr"];

  std::vector<unsigned int> numAtomsVec;
  int numSecondaryScatterers = jsonInput["nss"];
  if (numSecondaryScatterers == 0){ numAtomsVec = {npr}; }
  else { numAtomsVec = {npr,jsonInput["mss"]}; }

  //std::vector<unsigned int> numPrimarySecondaryAtoms {(unsigned int)(jsonInput["nss"]) };
  //if (numSeconaryScatterers > 0){ numPrimarySecondaryAtoms.push_back(unsigned int jsonInput["mss"])};

  //std::cout << "~~~~~~~  "<<numPrimarySecondaryAtoms[0] << std::endl;


  int    lat = jsonInput["lat"];
  int ilog = jsonInput["ilog"];
  int isym = 0;
  if (int(jsonInput["ncold"]) != 0){ isym = 1; }
  if (int(jsonInput["isabt"]) == 1){ isym += 2; }


  int iel = jsonInput["iel"];
  unsigned int secondaryScatterType = 0;
  if (numSecondaryScatterers > 0){ secondaryScatterType = jsonInput["b7"]; }

  double awr = jsonInput["awr"];
  int za = jsonInput["za"];
  double spr = jsonInput["spr"],
         sps = (numSecondaryScatterers == 0) ? 0.0 : double(jsonInput["sps"]);
  std::vector<double> temps;
  for (const auto& tempChunk : jsonInput["temperatures"]){
    temps.push_back(tempChunk["temperature"]);
  }

  double translationalWeight = jsonInput["temperatures"][temps.size()-1]["twt"];

  double aws = 0.0;
  if (numSecondaryScatterers > 0){ aws = jsonInput["aws"]; }
  std::vector<double> awrVec {awr,aws};
  if (numSecondaryScatterers == 0){ awrVec.resize(1); }


  //Float awr        = awrVec[0];
  //unsigned int npr = numAtomsVec[0];
  double sigma_b    = spr*pow(((1.0+awr)/awr),2);
  Range xsVec      = { spr*npr, sps };
  if (numSecondaryScatterers == 0){ xsVec.resize(1); }

  if (numSecondaryScatterers != 0 and secondaryScatterType <= 0){
    //Float aws = awrVec[1];
    double sigma_b2 = (aws == 0) ? 0 : sps*pow((1.0+aws)/aws,2);
    double srat=sigma_b2/sigma_b;
    for (size_t t = 0; t < temps.size(); ++t){
      for ( size_t a = 0; a < alphas.size(); ++a ){
        for ( size_t b = 0; b < betas.size(); ++b ){      
          sab[t][b+a*betas.size()] *= srat;
          sab[t][b+a*betas.size()] += principalScatterSAB[t][b+a*betas.size()];
        }
      }
    }
  }

  scaleDebyeWallerCoefficients( numSecondaryScatterers, secondaryScatterType, 
                                dwpix, dwp1, temps, awrVec );

  // write out the inelastic part
  auto epsilon = betas[betas.size()-1];
  auto emax    = 0.0253 * epsilon;
  unsigned int lasym = (isym > 1) ? 1 : 0;
  std::vector<unsigned int> secondaryScattererTypes {secondaryScatterType};
  if (numSecondaryScatterers == 0){ secondaryScattererTypes = {}; }

  ScatteringLawConstants constants(ilog, numSecondaryScatterers, epsilon, emax, 
    std::move(xsVec), std::move(awrVec), std::move(numAtomsVec), 
    std::move(secondaryScattererTypes));

  Inelastic mt4 = writeInelasticToENDF(sab,alphas,betas,temps,za,primaryTempf,
                                 secondaryTempf,lasym,lat,isym,ilog,constants);
  if (iel == 0 and translationalWeight == 0.0){
    // Write incoherent elastic part
    Elastic mt2(za,awr, writeIncElasticToENDF(sigma_b,temps,dwpix)); 
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }
  else if (iel > 0){
    // Write coherent elastic part
    Elastic mt2(za, awr, writeCohElasticToENDF( bragg, dwpix, dwp1, 
      numSecondaryScatterers, secondaryScatterType, numEdges, temps ));
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }

  return njoy::ENDFtk::file::Type<7>(std::move(mt4));

}











template <typename Range, typename Float>
auto endout( std::vector<Range>& sab, int za, Range awrVec, 
  const Float& spr, const Float& sps, const Range& temps, 
  int numSecondaryScatterers, unsigned int secondaryScatterType, 
  const std::vector<Range>& principalScatterSAB, const Range& alphas, 
  const Range& betas, Range& dwpix, Range& dwp1, int iel,
  const Float& translationalWeight, const Range& bragg, int numEdges, 
  Range primaryTempf, Range secondaryTempf, int ilog, int isym, int lat, 
  std::vector<unsigned int> numAtomsVec ){
  using std::pow;
  Float awr        = awrVec[0];
  unsigned int npr = numAtomsVec[0];
  Float sigma_b    = spr*pow(((1.0+awr)/awr),2);
  Range xsVec      = { spr*npr, sps };
  if (numSecondaryScatterers == 0){ xsVec.resize(1); }

  if (numSecondaryScatterers != 0 and secondaryScatterType <= 0){
    Float aws = awrVec[1];
    Float sigma_b2 = (aws == 0) ? 0 : sps*pow((1.0+aws)/aws,2);
    Float srat=sigma_b2/sigma_b;
    for (size_t t = 0; t < temps.size(); ++t){
      for ( size_t a = 0; a < alphas.size(); ++a ){
        for ( size_t b = 0; b < betas.size(); ++b ){      
          sab[t][b+a*betas.size()] *= srat;
          sab[t][b+a*betas.size()] += principalScatterSAB[t][b+a*betas.size()];
        }
      }
    }
  }

  scaleDebyeWallerCoefficients( numSecondaryScatterers, secondaryScatterType, 
                                dwpix, dwp1, temps, awrVec );

  // write out the inelastic part
  auto epsilon = betas[betas.size()-1];
  auto emax    = 0.0253 * epsilon;
  unsigned int lasym = (isym > 1) ? 1 : 0;
  std::vector<unsigned int> secondaryScattererTypes {secondaryScatterType};
  if (numSecondaryScatterers == 0){ secondaryScattererTypes = {}; }

  ScatteringLawConstants constants(ilog, numSecondaryScatterers, epsilon, emax, 
    std::move(xsVec), std::move(awrVec), std::move(numAtomsVec), 
    std::move(secondaryScattererTypes));

  Inelastic mt4 = writeInelasticToENDF(sab,alphas,betas,temps,za,primaryTempf,
                                 secondaryTempf,lasym,lat,isym,ilog,constants);
  if (iel == 0 and translationalWeight == 0.0){
    // Write incoherent elastic part
    Elastic mt2(za,awr, writeIncElasticToENDF(sigma_b,temps,dwpix)); 
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }
  else if (iel > 0){
    // Write coherent elastic part
    Elastic mt2(za, awr, writeCohElasticToENDF( bragg, dwpix, dwp1, 
      numSecondaryScatterers, secondaryScatterType, numEdges, temps ));
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }

  return njoy::ENDFtk::file::Type<7>(std::move(mt4));

}
