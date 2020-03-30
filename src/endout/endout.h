#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"
#include "coher/coher.h"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "endout/endout_util/writeIncElasticToENDF.h"
#include "endout/endout_util/writeCohElasticToENDF.h"
#include "ENDFtk.hpp"


using namespace njoy::ENDFtk;
using Elastic         = section::Type<7,2>;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using Inelastic       = section::Type<7,4>;
using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;


template <typename Range> 
auto scaleDebyeWallerCoefficients( int numSecondaryScatterers, 
  int secondaryScatterType, Range& dwpix, Range& dwp1, const Range& temps, 
  const Range& awrVec ){
  //const Float& awr, const Float& aws ){
  // display endf t-effective and debye-waller integral
  auto awr = awrVec[0];
  auto aws = awrVec[1];
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
auto endout( std::vector<Range>& sab, int za, Range awrVec, 
  const Float& spr, const Float& sps, const Range& temps, 
  int numSecondaryScatterers, unsigned int secondaryScatterType, 
  const std::vector<Range>& principalScatterSAB, const Range& alphas, 
  const Range& betas, Range& dwpix, Range& dwp1, int iel,
  const Float& translationalWeight, const Range& bragg, int numEdges, 
  Range tempf1, Range tempf, int ilog, int isym, int lat, 
  std::vector<unsigned int> numAtomsVec ){
  // compute bound scattering cross sections
  using std::pow;

  
  Float awr = awrVec[0];
  Float sigma_b  = spr*pow(((1.0+awr)/awr),2);
  Range xsVec = { spr, sps };
  if (numSecondaryScatterers == 0){ xsVec.resize(1); }
  Range primaryTempf   = (numSecondaryScatterers == 0) ? tempf1 : tempf ;
  Range secondaryTempf = (numSecondaryScatterers == 0) ? tempf  : tempf1;

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
