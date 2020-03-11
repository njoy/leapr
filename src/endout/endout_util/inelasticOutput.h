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


template <typename Float, typename Range, typename RangeOfRange, 
          typename ScatteringLawConstants >
auto writeToENDF( const RangeOfRange& fullSAB, const Range alphas, 
  const Range& betas, const Range& temps, const Float& za, Range effectiveTemps,
  int lasym, int lat, int isym, int ilog, ScatteringLawConstants constants ){

  using namespace njoy::ENDFtk;
  using ScatteringFunction     = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
  using ScatteringLaw          = section::Type< 7, 4 >::ScatteringLaw;
  using EffectiveTemperature   = section::Type< 7, 4 >::EffectiveTemperature;
  using Tabulated = section::Type< 7, 4 >::Tabulated;

  std::vector<Range> alphaVec (betas.size(),alphas);
  std::vector< long > boundaries   = { int(alphas.size()) },
                      interpolants = { int(temps.size() ) },
                      li (temps.size()-1,int(alphas.size()-1));

  std::vector<ScatteringFunction> chunkVectors;
  for (size_t b = 0; b < betas.size(); ++b){
    auto boundaries_b   = boundaries;
    auto interpolants_b = interpolants;
    auto li_b           = li;
    auto temps_b        = temps;
    auto alphas_b       = alphas;
    auto out_b          = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,b);
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

  EffectiveTemperature principal( { int(temps.size()) }, { int(temps.size()-1) }, 
                                  std::move(temps_c),
                                  std::move(effectiveTemps)
                                );

  section::Type< 7, 4 > chunk( za, awr, lat, lasym,
                               std::move( constants ),
                               std::move( scatter_law ),
                               std::move( principal ) );
  std::string buffer;
  auto output = std::back_inserter(buffer);
  chunk.print(output,27,7);
  std::cout << buffer << std::endl;

}




template <typename Range, typename RangeOfRange>
auto inelasticOutput( const Range& alphas, const Range& betas, const RangeOfRange& fullSAB,
  const Range& temps, int isym, int ilog, int lat ){
  using namespace njoy::ENDFtk;

  using namespace njoy::ENDFtk;
  //using ScatteringFunction     = section::Type< 7, 4 >::Tabulated::ScatteringFunction;


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

