#include <range/v3/all.hpp>
#include "generalTools/constants.h"
#include "generalTools/sigfig.h"
#include "ENDFtk.hpp"
#include <iostream>

template<typename Range, typename Float>
auto getlog10SAB( const Range& sab, const Float& be, int a, int b, int nbeta, int expBetaSign ){
  Float sabVal = sab[b+a*nbeta];
  if (sabVal <= 0.0){ return -999.0; }
  return sigfig(log(sabVal) + expBetaSign * be * 0.5, 7, 0);;
}

template<typename Range, typename Float>
auto getSAB( const Range& sab, const Float& be, int a, int b, int nbeta, int expBetaSign ){
  Float sabVal = sab[b+a*nbeta]*exp(expBetaSign*be*0.5);
  return (sabVal < 1e-9) ? sigfig(sabVal,6,0)
                         : sigfig(sabVal,7,0);
}

template <typename Range, typename RangeOfRange >
auto getSABreadyToWrite( const RangeOfRange& fullSAB, const Range& temps, 
  const Range& alphas, const Range& betas, int isym, int ilog, int lat, size_t b, 
  const RangeOfRange& fullSAB_2 = {std::vector<double> (0)} ){
    std::vector<Range> toWrite (temps.size()); 
    // This is a vector of vectors where the ith entry of this is a vector of 
    // SAB (in order of increasing temperature) for the ith temperature

    auto outputBeta = betas[b];
    size_t nbeta = betas.size();
    for (size_t t = 0; t < temps.size(); ++t){

      double sc = 1.0;
      if (lat == 1) {sc = 0.0253/(kb*temps[t]); }
      Range scr(alphas.size(),0.0);
      if ( isym == 0 or isym == 2){ outputBeta =  betas[b];         }
      else if ( b < nbeta )       { outputBeta = -betas[nbeta-b-1]; }
      else                        { outputBeta =  betas[b-nbeta+1]; }
      auto be = outputBeta*sc;
      
      for ( size_t a = 0; a < alphas.size(); ++a ){

        if (isym == 0 or isym == 2){
          int expBetaSign = ( isym == 0 ) ? -1 : 0;
          if (ilog == 0){ scr[a] = getSAB     (fullSAB[t],be,a,b,nbeta,expBetaSign);}
          else          { scr[a] = getlog10SAB(fullSAB[t],be,a,b,nbeta,expBetaSign);}
        } 

        if (isym == 1 or isym == 3){
          const Range& SAB = (b < nbeta-1) ? fullSAB[t] : fullSAB_2[t];
          int bIndex = std::abs(int(nbeta-b-1));
          int expBetaSign = ( isym == 1 ) ? 1 : 0; 
          if (ilog == 0){ scr[a] = getSAB     (SAB,be,a,bIndex,nbeta,expBetaSign);} 
          else          { scr[a] = getlog10SAB(SAB,be,a,bIndex,nbeta,expBetaSign);}
        }

        if (ilog == 0 and scr[a] < -999.0){ scr[a] = 0.0; }

      } // alpha loop
      toWrite[t] = scr;
    } // temperature loop

  return std::make_tuple(outputBeta,toWrite); 
}


template </*typename Float, typename ScatteringLawConstants,*/ typename Range, typename RangeOfRange >
auto writeToENDF( /*const Float& za, const Float& awr, int lasym, int lat, 
  const ScatteringLawConstants& constants,*/
  const RangeOfRange& fullSAB, const Range alphas, const Range& betas, const Range& temps 
  
  ){

  using namespace njoy::ENDFtk;
  using ScatteringFunction     = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
  using ScatteringLaw          = section::Type< 7, 4 >::ScatteringLaw;
  using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
  using EffectiveTemperature   = section::Type< 7, 4 >::EffectiveTemperature;
  //using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
  using Tabulated = section::Type< 7, 4 >::Tabulated;

  double za = 127.;
  double awr = 8.934780e+0;
  int lasym = 0;
  int lat = 1;
  ScatteringLawConstants constants( 0, 1.976285e+2, 5.000001e+0,
                                    6.153875e+0, 8.934780e+0, 1 );

  std::vector<Range> alphaVec (betas.size(),alphas);

  int isym = 0, ilog = 0;
  auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
  auto toWrite0 = std::get<1>(out0);

  auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
  auto toWrite1 = std::get<1>(out1);

  int nalpha = alphas.size();
  auto temps2  = temps;
  auto temps3  = temps;
  auto alphas2 = alphas;
  auto alphas3 = alphas;

  std::vector< long > boundaries = { 3 };
  std::vector< long > interpolants = { 2 };
  std::vector< long > li = { 4 };

  auto boundaries0 = boundaries;
  auto boundaries1 = boundaries;
  auto interpolants0 = interpolants;
  auto interpolants1 = interpolants;
  auto li0 = li;
  auto li1 = li;


  //std::vector<ScatteringFunction> chunkVectors(betas.size());


  ScatteringFunction chunk1 ( 
                              betas[0], 
                              std::move(boundaries0),
                              std::move(interpolants0),
                              std::move(temps2),
                              std::move(li0),
                              std::move(alphaVec[0]),
                              { std::move(toWrite0[0]),
                                std::move(toWrite0[1])
                              }
                              );
  ScatteringFunction chunk2 ( 
                              betas[1],
                              std::move(boundaries1),
                              std::move(interpolants1),
                              std::move(temps3),
                              std::move(li1),
                              std::move(alphaVec[1]),
                              { std::move(toWrite1[0]),
                                std::move(toWrite1[1])
                              }
                              );

  ScatteringLaw law =
    Tabulated( { 2 }, { 4 },
               { chunk1, 
                 chunk2
               } 
             );
  EffectiveTemperature principal( { 3 }, { 2 }, 
                                  { 293.6, 600., 1200. }, 
                                  { 5.332083e+2, 7.354726e+2,
                                    1.270678e+3 } );

      /*
  ScatteringLaw law =
    Tabulated( { 2 }, { 4 }, { 
      ScatteringFunction( 0.0, { 5 }, { 4 },                                       // beta0, 
              {296.0, 400.0},
        { 4 },
        { 0.1, 0.2, 0.3},
        { 
          { 2.386876e-4, 2.508466e-4, 2.636238e-4, 1.306574e-9, 5.29573e-10 },
          { 4.430020e-4, 4.655671e-4, 4.892796e-4, 4.510209e-8, 2.183942e-8 } 
        } ),
      ScatteringFunction( 0.2, { 5 }, { 2 },
              {296.0, 400.0},
        { 4 },
        { 0.1, 0.2, 0.3},
        { 
          { 2.386694e-4, 2.508273e-4, 2.636238e-4, 2.770291e-4, 2.911373e-4 },
          { 6.921141e-4, 7.273641e-4, 7.644060e-4, 8.033305e-4, 8.442328e-4 }
        } ) 
    } );

  EffectiveTemperature principal( { 3 }, { 2 }, 
                                  { 293.6, 600., 1200. }, 
                                  { 5.332083e+2, 7.354726e+2, 1.270678e+3 } );

                                  */

  /*
  ScatteringLaw law(
    Tabulated( { 2 }, { 4 },
      { ScatteringFunction(
          293.6, 0.0, { 5 }, { 4 },                                              // beta 1
          { 4.423802e-3, 4.649528e-3, 4.886772e-3, 8.418068e+1, 8.847604e+1 },   // alpha  vec 1
          { 2.386876e-4, 2.508466e-4, 2.636238e-4, 1.306574e-9, 5.29573e-10 } ), // sab for b1 alpha all
        ScatteringFunction(
          293.6, 3.952570e-2, { 5 }, { 2 },                                      // beta 2
          { 4.423802e-3, 4.649528e-3, 4.886772e-3, 8.418068e+1, 8.847604e+1 },   // alpha vec 2
          { 2.386694e-4, 2.508273e-4, 2.636238e-4, 2.770291e-4, 2.911373e-4 } )  // sab for b2 alpha all 
      } ) );
  EffectiveTemperature principal( { 3 }, { 2 }, { 293.6, 600., 1200. },
    { 5.332083e+2, 7.354726e+2, 1.270678e+3 } );
      */

  section::Type< 7, 4 > chunk( za, awr, lat, lasym,
                               std::move( constants ),
                               std::move( law ),
                               std::move( principal ) );
  std::string buffer;
  auto output = std::back_inserter(buffer);
  chunk.print(output,27,7);
  std::cout << buffer << std::endl;
  return;
  std::cout << fullSAB.size() << std::endl;
  std::cout << alphas.size() << std::endl;
  std::cout << betas.size() << std::endl;
  std::cout << temps.size() << std::endl;

}




template <typename Range, typename RangeOfRange>
auto inelasticOutput( const Range& alphas, const Range& betas, const RangeOfRange& fullSAB,
  const Range& temps, int isym, int ilog, int lat ){
  using namespace njoy::ENDFtk;

  using namespace njoy::ENDFtk;
  using ScatteringFunction     = section::Type< 7, 4 >::Tabulated::ScatteringFunction;


  //using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;
  //ScatteringLawConstants constants ( 0, 1, 1.9e2, 5.0, {6.15, 3.74}, {8.9, 1.5}, {1,2}, {1} );

  std::vector< long > boundaries = { 5 };
  std::vector< long > interpolants = { 4 };
  std::vector< long > li = { 4 };
  //std::vector< double > alphas =
  //{ 4.423802e-3, 4.649528e-3, 4.886772e-3, 8.418068e+1, 8.847604e+1 };
  std::vector< double > temperatures = { 293.6, 400.0 };
  std::vector< std::vector< double > > sab =
  { { 2.386876e-4, 2.508466e-4, 2.636238e-4, 1.306574e-9, 5.29573e-10 },
    { 4.430020e-4, 4.655671e-4, 4.892796e-4, 4.510209e-8, 2.183942e-8 } };
  //ScatteringFunction chunk( beta,
  //                          std::move( boundaries ),
  //                          std::move( interpolants ),
  //                          std::move( temperatures ),
  //                          std::move( li ),
  //                          std::move( alphas ),
  //                          std::move( sab ) );

  int nbt = betas.size();
  if (isym == 1 or isym == 3){ nbt = 2*betas.size()-1; }

  Range outputBetas (nbt,0.0);
  for (size_t b = 0; b < (unsigned) nbt; ++b){
    auto out = getSABreadyToWrite( fullSAB, temps, alphas, betas, isym, ilog, lat, b );
    outputBetas[b] = std::get<0>(out);
    auto toWrite   = std::get<1>(out);
    //for ( auto& vec : toWrite ){
    //  std::cout << (vec|ranges::view::all) << std::endl;
   // }
    //ScatteringFunction chunk(
    //                    betas[0], 
    //                    std::move(boundaries),
    //                    std::move(interpolants),
    //                    std::move(temperatures),
    //                    std::move(li),
    //                    std::move(alphas),
    //                    std::move(toWrite) );
                        //{ 3 }, 
                        //{ 4 },
                        //{ 293.6, 400 },
                        //{ 4 },
                        //std::copy(alphas),
                        //{ { 2.386876e-4, 2.508466e-4, 2.636238e-4 },
                        //  { 4.430020e-4, 4.655671e-4, 4.892796e-4 }
                        //} ),

    std::cout << std::endl;
    

  }
  
}

