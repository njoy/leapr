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







template <typename Range, typename RangeOfRange>
auto inelasticOutput( const Range& alphas, const Range& betas, const RangeOfRange& fullSAB,
  const Range& temps, int isym, int ilog, int lat ){
  using namespace njoy::ENDFtk;


  using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;
  ScatteringLawConstants constants ( 0, 1, 1.9e2, 5.0, {6.15, 3.74}, {8.9, 1.5}, {1,2}, {1} );


  int nbt = betas.size();
  if (isym == 1 or isym == 3){ nbt = 2*betas.size()-1; }

  Range outputBetas (nbt,0.0);
  for (size_t b = 0; b < (unsigned) nbt; ++b){
    auto out = getSABreadyToWrite( fullSAB, temps, alphas, betas, isym, ilog, lat, b );
    outputBetas[b] = std::get<0>(out);
    auto toWrite = std::get<1>(out);
    for ( auto& vec : toWrite ){
        std::cout << (vec|ranges::view::all) << std::endl;
    }
    std::cout << std::endl;
    

  }
  
}

