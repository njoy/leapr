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
                 std::ostream& ,
                 const nlohmann::json& ){

  //output << "Input arguments:\n" << jsonInput.dump(2) << std::endl;

  std::setprecision(15);
  int numSecondaryScatterers = jsonInput["nss"];
  int b7 = 0;
  if (numSecondaryScatterers > 0){ b7 = jsonInput["b7"]; }
  int numIter = 2;
  if ( numSecondaryScatterers == 0 or b7 > 0 ){ numIter = 1; }


  double awr = jsonInput["awr"];
  double aws = 0.0;
  if (numSecondaryScatterers > 0){ aws = jsonInput["aws"]; }

  auto temperatureInfo = jsonInput["temperatures"];
  int ntempr = jsonInput["ntempr"];
  std::vector<double> dwpix (ntempr,0.0),
                      dwp1  (ntempr,0.0),
                      primaryEffectiveTemp(ntempr,0.0),
                      secondaryEffectiveTemp(ntempr,0.0),
                      primaryDWF(ntempr,0.0),
                      secondaryDWF(ntempr,0.0);

  int ncold = jsonInput["ncold"];

  const std::vector<double> alphas = jsonInput["alpha"],
                            betas  = jsonInput["beta"];

  std::vector<std::vector<double>> primarySabAllTemps(ntempr);
  std::vector<std::vector<double>> secondarySabAllTemps(ntempr);
  std::vector<std::vector<double>> primarySabAllTempsPosBeta;
  if (ncold > 0){ primarySabAllTempsPosBeta.resize(ntempr); }

  std::vector<double> temperatures(ntempr);

  for (int scatterIter = 0; scatterIter < numIter; ++scatterIter){

    int secScatterOffset = (scatterIter == 1 and b7 == 0) ?
        int(jsonInput["ntempr"]) : 0;

    for (int itemp = 0; itemp < int(jsonInput["ntempr"]); ++itemp){

      std::vector<double> sab(alphas.size()*betas.size(),0.0);
      std::vector<double> sab2(alphas.size()*betas.size(),0.0);

      auto tempInfo = temperatureInfo[itemp+secScatterOffset];
      double temperature = tempInfo["temperature"];
      output << std::fixed << std::setprecision(2) << "\n\nBeginning Calculation for " << 
                temperature << " K " << std::endl << std::endl;
      if (scatterIter == 0){ temperatures[itemp] = temperature; }

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
                           continuousWgt, rho, scaledAlphas, scaledBetas, sab,
                           output, int(jsonInput["iprint"]) );

      dwpix[itemp] = std::get<0>(continOutput);
      double effectiveTemperature = std::get<1>(continOutput)*temperature;
      output << "            - Debye-Waller Factor       " << dwpix[itemp]<< std::endl << 
                "            - Effective Temperature     " << effectiveTemperature<< std::endl;     


      double transWgt = tempInfo["twt"];
      if (transWgt > 0){
        output << "\n\n     Beginning translational calculation... " << std::endl;
        translational( scaledAlphas, scaledBetas, transWgt, rho_dx,
                       static_cast<double>(tempInfo["c"]), dwpix[itemp], continuousWgt,
                       effectiveTemperature , temperature, sab );
      }

      if (not tempInfo["oscillators"]["energies"].is_null()){
        std::vector<double> oscE = tempInfo["oscillators"]["energies"],
                            oscW = tempInfo["oscillators"]["weights"];
        output << "\n\n     Beginning discrete oscillator calculation... " << std::endl;
        discreteOscillators( dwpix[itemp], transWgt, continuousWgt, scaledAlphas,
            scaledBetas, temperature, ranges::view::zip(oscE,oscW),
            effectiveTemperature, sab );
      }

      if (ncold > 0){
        bool free = false;
        auto pairCorrelationInfo = tempInfo["pairCorrelation"];
        std::vector<double> kappa = pairCorrelationInfo["skappa"];
        double                dka = pairCorrelationInfo["dka"];
        double tbart = effectiveTemperature/temperature;
        output << "\n\n     Beginning cold hydrogen calculation... " << std::endl;
        coldHydrogen( tev, ncold, transWgt+continuousWgt, alphas, betas, dka,
                      kappa, free, sab, sab2, tbart );
        primarySabAllTempsPosBeta[itemp] = sab2;
      }

      if (int(jsonInput["nsk"]) == 2 and ncold == 0){
        std::vector<double> skappa = tempInfo["pairCorrelation"]["skappa"];
        double dka = tempInfo["pairCorrelation"]["dka"];
        double cfrac = tempInfo["pairCorrelation"]["cfrac"];
        skold( cfrac, tev, scaledAlphas, betas, skappa, awr, dka, sab );
      } 


      if (scatterIter == 1 and b7 == 0){
        secondaryEffectiveTemp[itemp] = effectiveTemperature ;
        secondarySabAllTemps[itemp] = sab;
        secondaryDWF[itemp] = dwpix[itemp];
      }
      else{
        primaryEffectiveTemp[itemp] = effectiveTemperature ;
        primarySabAllTemps[itemp] = sab;
        primaryDWF[itemp] = dwpix[itemp];
      }

    } // temp loop

  } // Primary and Secondary Scatter Loop

  //---------------- Coherent (Elastic) ----------------------
  int iel = jsonInput["iel"];
  unsigned int npr = jsonInput["npr"];
  std::vector<double> bragg ( 60000, 0.0 );
  int nedge = 0;
  if (iel > 0){
    double emax = 5.0;
    output << "\n\n     Beginning coherent elastic calculation... " << std::endl;
    auto coherentElasticOut = coherentElastic( iel, npr, bragg, emax );
    nedge = std::get<1>(coherentElasticOut)* 0.5;
    bragg.resize(std::get<0>(coherentElasticOut));
  }
  else {
    bragg.resize(0);
  }


  //---------------------------- Write Output --------------------------------
  output << "\n\n     Writing thermal scattering information to tape" << std::endl;
  writeOutput( jsonInput, primarySabAllTemps, secondarySabAllTemps, 
   primaryDWF, secondaryDWF, bragg, nedge,  primaryEffectiveTemp,  
   secondaryEffectiveTemp, primarySabAllTempsPosBeta, temperatures );



}
};

}
}
