#include "continuous/continuous.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include "coher/coher.h"
#include "skold/skold.h"
#include "generalTools/constants.h"
#include <range/v3/all.hpp>
#include <variant>
#include "endout/endout.h"
//#include "lipservice.hpp"


template <typename Range, typename RangeOfRanges, typename Float>
auto leapr( std::vector<int> generalInfo, std::vector<int> scatterControl,
  Range scatterInfo, Range temps, Range alphas, Range betas, RangeOfRanges rhoVec, 
  Range rho_dx_vec, RangeOfRanges transInfo, RangeOfRanges oscE_vec, 
  RangeOfRanges oscW_vec, Float smin, std::tuple<std::vector<double>,double> 
  kappaInfo = std::tuple<std::vector<double>,double>() ){

  std::cout.precision(15);
  int nphon = generalInfo[0];

  // Do we have a secondary scatterer?
  int numSecondaryScatterers = scatterControl[4],
                          b7 = scatterControl[5];
  int numIter = 2;
  if ( numSecondaryScatterers == 0 or b7 > 0 ){ numIter = 1; }


  int  lat = generalInfo[5];
  auto awr = scatterInfo[0];
  auto aws = scatterInfo[2];

  Range dwpix (temps.size(),0.0),
        tempf (temps.size(),0.0),
        dwp1  (temps.size(),0.0),
        tempf1(temps.size(),0.0);

  int ncold = scatterControl[2];
  int isabt = generalInfo[3];

  std::vector<std::vector<Float>> sab_AllTemps(temps.size());
  for (int scatterIter = 0; scatterIter < numIter; ++scatterIter){

    for (size_t itemp = 0; itemp < temps.size(); ++itemp){

      std::vector<Float> sab(alphas.size()*betas.size(),0.0);
      std::vector<Float> sab2(alphas.size()*betas.size(),0.0);
      auto tev = kb * temps[itemp];

      Float sc   = (lat         == 1) ? 0.0253/tev : 1.0,
            arat = (scatterIter == 0) ?        1.0 : aws/awr;
      Float scaling = sc/arat;
      Range scaledAlphas = alphas,
            scaledBetas  = betas;
      for ( auto& a : scaledAlphas ){ a *= scaling; }
      for ( auto& b : scaledBetas  ){ b *= sc;      }


      //---------------- Incoherent (Elastic and Inelastic) ----------------------
      auto rho_dx    = rho_dx_vec[itemp]/tev;
      auto continWgt = transInfo[2][itemp];
      auto rho       = rhoVec[itemp];
      auto continOutput = contin(nphon, rho_dx, continWgt, rho, scaledAlphas, scaledBetas, sab);

      dwpix[itemp] = std::get<0>(continOutput);
      tempf[itemp] = std::get<1>(continOutput)*temps[itemp];
    
      auto transWgt = transInfo[0][itemp];
      if (transWgt > 0){
        auto diffusion_const = transInfo[1][itemp];
        trans( scaledAlphas, scaledBetas, transWgt, rho_dx, diffusion_const, 
               dwpix[itemp], continWgt, tempf[itemp], temps[itemp], sab );
      }

      auto oscEnergiesWgts = ranges::view::zip(oscE_vec[itemp],oscW_vec[itemp]);
      if (oscEnergiesWgts.size() > 0){
        discre( dwpix[itemp], transWgt, continWgt, scaledAlphas, scaledBetas, 
                temps[itemp], oscEnergiesWgts, tempf[itemp], sab );
      }

      if (ncold > 0){
        bool free = false;

        std::vector<double> kappa = std::get<0>(kappaInfo);
        double dka                = std::get<1>(kappaInfo);
        coldh( tev, ncold, transWgt+continWgt, alphas, betas, dka, kappa, free, 
               sab, sab2, tempf[itemp]);
        //std::cout << (sab|ranges::view::all) << std::endl;
        //std::cout << (sab2|ranges::view::all) << std::endl;
      }


      sab_AllTemps[itemp] = sab;

    } // temp loop

    if (scatterIter == 0 and numIter == 2){
      for (size_t itemp = 0; itemp < temps.size(); ++itemp){
        tempf1[itemp] = tempf[itemp];
        dwp1[itemp]   = dwpix[itemp];
      }
    }
    //std::cout << (dwp1|ranges::view::all) << std::endl;
    //std::cout << (tempf1|ranges::view::all) << std::endl;

  } // Primary and Secondary Scatter Loop
  

  //---------------- Coherent (Elastic) ----------------------
  int iel = scatterControl[1];
  unsigned int npr = scatterControl[0];
  std::variant<Range,bool> braggOutput;
  std::vector<double> bragg ( 60000, 0.0 );
  int nedge = 0;
  if (iel > 0){
    double emax = 5.0;
    auto coherOut = coher( iel, npr, bragg, emax );
    //int numVals = std::get<0>(coherOut);
    //int numVals2 = std::get<1>(coherOut);
    nedge = std::get<1>(coherOut)* 0.5;
    //std::cout << numVals << "   " << numVals2 << "   " << nedge  << std::endl;
    bragg.resize(std::get<0>(coherOut)); 
    braggOutput = bragg;
  }
  else { 
      braggOutput = false; 
      bragg.resize(0);
  }
  int numEdges = bragg.size();
  //std::cout << numEdges << std::endl;

  //try {
  //  std::get<bool>(braggOutput); // w contains int, not float: will throw
  //  //std::cout << "No bragg!" << std::endl;
  //}
  //catch (const std::bad_variant_access&) {}



  //---------------------------- Write Output --------------------------------

  int za = generalInfo[2];
  //std::vector<Float> awrVec = (numSecondaryScatterers == 0) ? {awr} : {awr,aws};
  std::vector<Float> awrVec {awr,aws};
  if (numSecondaryScatterers == 0){ awrVec.resize(1); }

  Float spr = scatterInfo[1],
        sps = scatterInfo[3];
  
  int ilog = generalInfo[4];


  int isym = 0;
  if (ncold != 0){ isym = 1; }
  if (isabt == 1){ isym += 2; }


  auto transWgt = transInfo[0][temps.size()-1];
  unsigned int mss = scatterControl[6];
  std::vector<unsigned int> numAtomsVec {npr,mss};
  if (numSecondaryScatterers == 0){ numAtomsVec.resize(1); }

    njoy::ENDFtk::file::Type<7> MF7 = endout(sab_AllTemps,za,awrVec,spr,sps,temps,
          numSecondaryScatterers, b7 , sab_AllTemps,alphas,
          betas,dwpix,dwp1,iel,transWgt,bragg,nedge,tempf,tempf1,ilog,
          isym,lat,numAtomsVec);

  return MF7;
  /*
return;
  */
  std::cout << numEdges << std::endl;
  std::cout << ilog<< std::endl;
  std::cout << za<< std::endl;
  /*
  */
  

  /*
  //njoy::ENDFtk::section::Type<7,2> MT2 = MF7.MT(2_c);
  //njoy::ENDFtk::section::Type<7,4> MT4 = MF7.MT(4_c);

  */
  //std::cout << MT4.ZA() << std::endl;

  //std::string buffer;
  //auto output = std::back_inserter( buffer );
  //MT4.print( output, 27, 7 );
  //std::cout << buffer << std::endl;


        //std::string buffer;
        //auto output = std::back_inserter( buffer );
        //chunk.print( output, 27, 7 );




  /*
    //return;
    std::cout << ilog << sps << spr << za << std::endl;
    std::cout << generalInfo.size() << std::endl;
    std::cout << scatterControl.size() << std::endl;
    std::cout << scatterInfo.size() << std::endl;
    std::cout << temps.size() << std::endl;
    std::cout << alphas.size() << std::endl;
    std::cout << betas.size() << std::endl;
    std::cout << rhoVec.size() << std::endl;
    std::cout << rho_dx_vec.size() << std::endl;
    std::cout << transInfo.size() << std::endl;
    std::cout << oscE_vec.size() << std::endl;
    std::cout << oscW_vec.size() << std::endl;
    */
    std::cout << smin << std::endl;

}


        /*
template<typename Range, typename Float, typename RangeZipped>
auto bigFunc( int nphon, Float awr, int iel, int npr, int ncold, Float aws, int lat, Range alpha, Range beta, Range temps, Float delta, Range rho, Float transWgt, Float diffusion_const, Float continWgt, RangeZipped oscEnergiesWgts, Float dka, Range kappaVals, Float cfrac, std::tuple<int,int,Range> secondaryScatterInput ){

  auto bigTuple = leapr(nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temps, delta, rho, transWgt, diffusion_const, continWgt, oscEnergiesWgts, dka, kappaVals, cfrac, secondaryScatterInput );

  auto sab             = std::get<0>(bigTuple);
  auto effectiveTemps  = std::get<1>(bigTuple);
  auto lambdaVals      = std::get<2>(bigTuple);
  auto braggOutput     = std::get<3>(bigTuple);
  auto sab2            = std::get<4>(bigTuple);
  auto effectiveTemps2 = std::get<5>(bigTuple);
  auto lambdaVals2     = std::get<6>(bigTuple);
  auto braggOutput2    = std::get<7>(bigTuple);
  auto sab_AllTemps    = std::get<8>(bigTuple);

  std::cout << (sab_AllTemps[0]|ranges::view::all) << std::endl;
  std::vector<Float> awrVec {awr,aws};


  //njoy::ENDFtk::file::Type<7> MF7 = endout(sab_AllTemps,za,awrVec,spr,sps,temps,nuumSecondaryScatterers,secondaryScattererType,principalScatterSAB,alphas,betas,dwpix,dwp1,iel,translationalWeight,bragg,numEdges,tempf,tempf1,ilog,isym,lat,numAtomsVec);
//          std::make_tuple(sab, effectiveTemps, lambdaVals, braggOutput,
//                         sab2,effectiveTemps2,lambdaVals2,braggOutput2);





}

*/








