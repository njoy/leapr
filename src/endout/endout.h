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
auto scaleDebyeWallerCoefficients( const nlohmann::json& jsonInput,
  Range& primaryDWF, Range& secondaryDWF, const Range& temps, 
  const Range& awrVec ){

  int numSecondaryScatterers = jsonInput["nss"];
  unsigned int secondaryScatterType = 0;
  if (numSecondaryScatterers > 0){ secondaryScatterType = jsonInput["b7"]; }

  auto awr = awrVec[0], aws = awrVec[1];

  for (size_t i = 0; i < temps.size(); ++i){
    if (numSecondaryScatterers == 0 or secondaryScatterType > 0){
       primaryDWF[i] /= (awr*temps[i]*kb);
    }
    else {
       primaryDWF[i] /= (aws*temps[i]*kb);
       secondaryDWF[i]  /= (awr*temps[i]*kb);
    }
  }
}


template <typename Range>
auto endout( const nlohmann::json& jsonInput, std::vector<Range>& sab,
  const std::vector<Range>& secondaryScatterSAB, const Range& alphas, 
  const Range& betas, Range& primaryDWF, Range& secondaryDWF, const Range& bragg, 
  int numEdges, Range primaryTempf, Range secondaryTempf, const Range& temps){ 

  using std::pow;

  std::vector<unsigned int> numAtomsVec;
  unsigned int npr = jsonInput["npr"];
  int numSecondaryScatterers = jsonInput["nss"];
  if (numSecondaryScatterers == 0){ numAtomsVec = {npr}; }
  else { numAtomsVec = {npr,jsonInput["mss"]}; }

  int isym = 0;
  if (int(jsonInput["ncold"]) != 0){ isym = 1; }
  if (int(jsonInput["isabt"]) == 1){ isym += 2; }


  int iel = jsonInput["iel"];
  unsigned int secondaryScatterType = 0;
  if (numSecondaryScatterers > 0){ secondaryScatterType = jsonInput["b7"]; }

  int za = jsonInput["za"];
  double spr = jsonInput["spr"],
         sps = (numSecondaryScatterers == 0) ? 0.0 : double(jsonInput["sps"]);

  int ntempr = jsonInput["ntempr"];

  double translationalWeight = jsonInput["temperatures"][temps.size()-1]["twt"];
  double awr = jsonInput["awr"];
  double aws = (numSecondaryScatterers > 0) ? double(jsonInput["aws"]) : 0.0;
  std::vector<double> awrVec {awr,aws};

  double sigma_b    = spr*pow(((1.0+awr)/awr),2);
  Range xsVec      = { spr*npr, sps };
  if (numSecondaryScatterers == 0){ xsVec.resize(1); }

  if (numSecondaryScatterers != 0 and secondaryScatterType <= 0){
    double sigma_b2 = (aws == 0) ? 0 : sps*pow((1.0+aws)/aws,2);
    double srat=sigma_b2/sigma_b;

    for (size_t t = 0; t < ntempr; ++t){
      for ( size_t a = 0; a < alphas.size(); ++a ){
        for ( size_t b = 0; b < betas.size(); ++b ){      
          sab[t][b+a*betas.size()] += srat*secondaryScatterSAB[t][b+a*betas.size()];
        }
      }
    }
  }

  scaleDebyeWallerCoefficients( jsonInput, primaryDWF, secondaryDWF, temps, awrVec );

  // write out the inelastic part
  auto epsilon = betas[betas.size()-1];
  auto emax    = 0.0253 * epsilon;
  unsigned int lasym = (isym > 1) ? 1 : 0;
  std::vector<unsigned int> secondaryScattererTypes {secondaryScatterType};
  if (numSecondaryScatterers == 0){ secondaryScattererTypes = {}; }
  if (numSecondaryScatterers == 0){ awrVec.resize(1); }
  int ilog = jsonInput["ilog"];
  ScatteringLawConstants constants(ilog, numSecondaryScatterers, epsilon, emax, 
    std::move(xsVec), std::move(awrVec), std::move(numAtomsVec), 
    std::move(secondaryScattererTypes));

  Inelastic mt4 = writeInelasticToENDF(sab,alphas,betas,temps,za,primaryTempf,
                                 secondaryTempf,lasym,int(jsonInput["lat"]),isym,ilog,constants);
  if (iel == 0 and translationalWeight == 0.0){
    // Write incoherent elastic part
    Elastic mt2(za,awr, writeIncElasticToENDF(sigma_b,temps,primaryDWF)); 
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }
  else if (iel > 0){
    // Write coherent elastic part
    Elastic mt2(za, awr, writeCohElasticToENDF( bragg, primaryDWF, secondaryDWF, 
      numSecondaryScatterers, secondaryScatterType, numEdges, temps ));
    return njoy::ENDFtk::file::Type<7>( std::move(mt2), std::move(mt4) );
  }

  return njoy::ENDFtk::file::Type<7>(std::move(mt4));

}








void writeOutput( const nlohmann::json jsonInput, 
  std::vector<std::vector<double>> primarySabAllTemps,
  std::vector<std::vector<double>> secondarySabAllTemps,
  std::vector<double> primaryDWF, std::vector<double> secondaryDWF,
  std::vector<double> bragg, int nedge, std::vector<double> tempf, std::vector<double> tempf1,
  std::vector<double> temperatures ){

  const std::vector<double> alphas = jsonInput["alpha"],
                            betas  = jsonInput["beta"];

  njoy::ENDFtk::file::Type<7> MF7 = endout( jsonInput, primarySabAllTemps, 
    secondarySabAllTemps, alphas, betas, primaryDWF, secondaryDWF, bragg, nedge, tempf, tempf1, temperatures );

  std::vector< njoy::ENDFtk::DirectoryRecord > index = {};
  if ( MF7.hasSection( 2 ) ) {
    index.emplace_back( 7, 2, MF7.section( 2_c ).NC(), 0 );
  }
  index.emplace_back( 7, 4, MF7.section( 4_c ).NC(), 0 );

  std::string buffer2;
  auto output2 = std::back_inserter( buffer2 );
  auto dir = index[0];
  dir.print( output2, 1, 1, 451 );


  int    lrp  = -1, lfi  = 0, nlib = 0, nmod = 0;
  double elis =  0, sta  = 0;
  int    lis  =  0, liso = 0, nfor = 6;
  double awi  =  1, emax = 0;
  int    lrel =  0, nsub = 0, nver = 0;
  double temp =  0;
  int    ldrv =  0;

  std::string comments;
  for (std::string comment : jsonInput["comments"]){
    comments = comments + comment;
  }
  double awr = jsonInput["awr"];

  double zaid = jsonInput["za"];
  njoy::ENDFtk::section::Type< 1, 451 > mf1mt451(
    zaid, awr, lrp, lfi, nlib, nmod, elis, sta, lis, liso, nfor, awi, emax, 
    lrel, nsub, nver, temp, ldrv, std::move(comments), std::move( index ) );


  njoy::ENDFtk::Material material( int(jsonInput["mat"]), 
                          njoy::ENDFtk::file::Type<1>( std::move( mf1mt451 ) ), 
                          std::move( MF7 ) );
  std::vector<njoy::ENDFtk::Material> materials {material};

  std::string title = jsonInput["title"];
  njoy::ENDFtk::Tape tape( TapeIdentification(std::move(title),1), std::move(materials) );

  std::string buffer;
  auto materialOutput = std::back_inserter( buffer );
  tape.print( materialOutput );
  int nout = jsonInput["nout"];
  std::string tapeNumber = std::to_string(nout);
  std::string name = "tape"+tapeNumber;
  std::ofstream out(name);
  out << buffer;
  out.close();

}






