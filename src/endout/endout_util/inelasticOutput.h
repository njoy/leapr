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

/*

template <typename Float, typename Range, typename RangeOfRange, 
          typename ScatteringLawConstants >
auto writeToENDF( const RangeOfRange& fullSAB, const Range alphas, 
  const Range& betas, const Range& temps, const Float& za, 
  Range effectiveTempsPrincipal, Range effectiveTempsSecondary,
  int lasym, int lat, int isym, int ilog, ScatteringLawConstants constants ){

  using namespace njoy::ENDFtk;
  using ScatteringFunction   = section::Type<7,4>::Tabulated::ScatteringFunction;
  using ScatteringLaw        = section::Type<7,4>::ScatteringLaw;
  using EffectiveTemperature = section::Type<7,4>::EffectiveTemperature;
  using Tabulated            = section::Type<7,4>::Tabulated;

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

  EffectiveTemperature principal( { int(temps.size()) }, { int(temps.size()-1) }, 
                                  std::move(temps_c),
                                  std::move(effectiveTempsPrincipal)
                                );
  EffectiveTemperature secondary( { int(temps.size()) }, { int(temps.size()-1) }, 
                                  std::move(temps_d),
                                  std::move(effectiveTempsSecondary)
                                );

  section::Type< 7, 4 > chunk( za, awr, lat, lasym,
                               std::move(constants),
                               std::move(scatter_law),
                               std::move(principal), 
                             { std::optional<EffectiveTemperature>(secondary) } );
  return chunk;
}
*/
