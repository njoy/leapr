#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"
#include "coher/coher.h"
#include "ENDFtk.hpp"


using namespace njoy::ENDFtk;
using Elastic         = section::Type<7,2>;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using Inelastic       = section::Type<7,4>;


template <typename Range, typename Float> 
auto scaleDebyeWallerCoefficients( int numSecondaryScatterers, 
  int secondaryScatterType, Range& dwpix, Range& dwp1, const Range& temps, 
  const Float& awr, const Float& aws ){
  // display endf t-effective and debye-waller integral
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





template <typename Range, typename Float>
auto endout( Range& sab, const Float awr, const Float& aws, const Float& spr, 
  const Float& sps, const Range& temps, int numSecondaryScatterers, 
  int secondaryScatterType, const Range& secondaryScatterVecThing, 
  const Range& alphas, const Range& betas, Range& dwpix, Range& dwp1, int iel,
  const Float& translationalWeight, const Range& bragg ){
  // compute bound scattering cross sections
  using std::pow;
  Float sb  = spr*pow(((1.0+awr)/awr),2);
  Float sbs = (aws == 0) ? 0 : sps*pow((1.0+aws)/aws,2);

  // for mixed moderators, merge ssm results
  if (numSecondaryScatterers != 0 and secondaryScatterType <= 0){
    Float srat=sbs/sb;
    for (size_t t = 0; t < temps.size(); ++t){
      for ( size_t a = 0; a < alphas.size(); ++a ){
        for ( size_t b = 0; b < betas.size(); ++b ){      
          sab[b+a*betas.size()] = srat*sab[b+a*betas.size()] 
                                + secondaryScatterVecThing[b+a*betas.size()];
        }
      }
    }
  }

  scaleDebyeWallerCoefficients( numSecondaryScatterers, secondaryScatterType, 
                                dwpix, dwp1, temps, awr, aws );


  if (iel == 0 and translationalWeight == 0.0){
    iel = -1;
    // Write out the incoherent elastic part
    //Elastic 
    

  }
  else if (iel > 0){
    // Write out the coherent elastic part
    //auto cohElasticChunk = writeCohElasticToENDF( bragg, dwpix, dwp1, 
    //  numSecondaryScatterers, secondaryScatterType, numEdges, tol, temps );
      std::cout << bragg.size() << std::endl;


  
  }


  // write out the inelastic part


  return;
  std::cout << sab[0] << std::endl;
}
