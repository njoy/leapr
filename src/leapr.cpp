#include "continuous/continuous.h"
#include "translational/translational.h"
#include "discreteOscillators/discreteOscillators.h"
#include "coldHydrogen/coldHydrogen.h"
#include "coherentElastic/coherentElastic.h"
#include "skold/skold.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>
#include <variant>
#include "endout/endout.h"
#include <string>
using namespace njoy;
#include "lipservice.hpp"

namespace njoy {
namespace LEAPR {

class LEAPR{

public:
void operator()( const nlohmann::json& jsonInput,
                 std::ostream& output,
                 std::ostream& error,
                 const nlohmann::json& ){

  output << "Input arguments:\n" << jsonInput.dump(2) << std::endl;
  output << fmt::format( "Input arguments:\n{}", jsonInput.dump(2) ) << std::endl;

  // Do we have a secondary scatterer?
  int numSecondaryScatterers = jsonInput["nss"];
  int b7 = 0;
  if (numSecondaryScatterers > 0){ b7 = jsonInput["b7"]; }
  int numIter = 2;
  if ( numSecondaryScatterers == 0 or b7 > 0 ){ numIter = 1; }


  double awr = jsonInput["awr"];
  double aws = 0.0;
  if (numSecondaryScatterers > 0){ aws = jsonInput["aws"]; }

  auto temperatureInfo = jsonInput["temperatures"];
  std::vector<double> dwpix (temperatureInfo.size(),0.0),
                      tempf (temperatureInfo.size(),0.0),
                      dwp1  (temperatureInfo.size(),0.0),
                      tempf1(temperatureInfo.size(),0.0);

  int ncold = jsonInput["ncold"];

  std::vector<double> alphas = jsonInput["alpha"],
                      betas  = jsonInput["beta"];

  std::vector<std::vector<double>> sab_AllTemps(temperatureInfo.size());
  std::vector<double> temperatures(temperatureInfo.size());
  for (int scatterIter = 0; scatterIter < numIter; ++scatterIter){

    for (size_t itemp = 0; itemp < temperatureInfo.size(); ++itemp){

      auto tempInfo = temperatureInfo[itemp];

      double temperature = tempInfo["temperature"];
      if (scatterIter == 0){ temperatures[itemp] = temperature; }

      std::vector<double> sab(alphas.size()*betas.size(),0.0);
      std::vector<double> sab2(alphas.size()*betas.size(),0.0);
      double tev = kb * temperature;

      int    lat  = jsonInput["lat"];
      double sc   = (lat         == 1) ? 0.0253/tev : 1.0,
             arat = (scatterIter == 0) ?        1.0 : aws/awr;
      double scaling = sc/arat;
      std::vector<double> scaledAlphas = alphas,
                          scaledBetas  = betas;
      for ( double& a : scaledAlphas ){ a *= scaling; }
      for ( double& b : scaledBetas  ){ b *= sc;      }


      //---------------- Incoherent (Elastic and Inelastic) ----------------------
      double rho_dx           = double(tempInfo["delta"])/tev;
      double continuousWgt    = tempInfo["tbeta"];
      std::vector<double> rho = tempInfo["rho"];
      auto continOutput  = continuous(int(jsonInput["nphon"]), rho_dx,
                           continuousWgt, rho, scaledAlphas, scaledBetas, sab);
      dwpix[itemp] = std::get<0>(continOutput);
      tempf[itemp] = std::get<1>(continOutput)*temperature;

      double transWgt = tempInfo["twt"];
      if (transWgt > 0){
        translational( scaledAlphas, scaledBetas, transWgt, rho_dx,
                       static_cast<double>(tempInfo["c"]), dwpix[itemp], continuousWgt,
                       tempf[itemp], temperature, sab );
      }

      if (not tempInfo["oscillators"]["energies"].is_null()){
        std::vector<double> oscE = tempInfo["oscillators"]["energies"],
                            oscW = tempInfo["oscillators"]["weights"];
        discreteOscillators( dwpix[itemp], transWgt, continuousWgt, scaledAlphas,
            scaledBetas, temperature, ranges::view::zip(oscE,oscW),
            tempf[itemp], sab );
      }

      if (ncold > 0){
        bool free = false;
        auto pairCorrelationInfo = tempInfo["pairCorrelation"];
        std::vector<double> kappa = pairCorrelationInfo["skappa"];
        double                dka = pairCorrelationInfo["dka"];
        coldHydrogen( tev, ncold, transWgt+continuousWgt, alphas, betas, dka,
                      kappa, free, sab, sab2, tempf[itemp]);
      }

      sab_AllTemps[itemp] = sab;

    } // temp loop

    if (scatterIter == 0 and numIter == 2){
      for (size_t itemp = 0; itemp < temperatureInfo.size(); ++itemp){
        tempf1[itemp] = tempf[itemp];
        dwp1[itemp]   = dwpix[itemp];
      }
    }

  } // Primary and Secondary Scatter Loop


  //---------------- Coherent (Elastic) ----------------------
  int iel = jsonInput["iel"];
  unsigned int npr = jsonInput["npr"];
  std::variant<std::vector<double>,bool> braggOutput;
  std::vector<double> bragg ( 60000, 0.0 );
  int nedge = 0;
  if (iel > 0){
    double emax = 5.0;
    auto coherentElasticOut = coherentElastic( iel, npr, bragg, emax );
    nedge = std::get<1>(coherentElasticOut)* 0.5;
    bragg.resize(std::get<0>(coherentElasticOut));
    braggOutput = bragg;
  }
  else {
    braggOutput = false;
    bragg.resize(0);
  }

  //---------------------------- Write Output --------------------------------

  njoy::ENDFtk::file::Type<7> MF7 = endout( jsonInput, sab_AllTemps, 
    sab_AllTemps, alphas, betas, dwpix, dwp1, bragg, nedge, tempf, tempf1 );

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

  double zaid = jsonInput["za"];
  njoy::ENDFtk::section::Type< 1, 451 > mf1mt451(
    zaid, awr, lrp, lfi, nlib, nmod, elis, sta, lis, liso, nfor, awi, emax, 
    lrel, nsub, nver, temp, ldrv, std::move(comments), std::move( index ) );

  //output << "this is output";
  //error << "this is error";

  int mat = jsonInput["mat"];

  njoy::ENDFtk::Material material( mat, 
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
};

}
}
