#include "catch.hpp"
#include "generalTools/testing.h"
#include "endout/endout_util/writeIncElasticToENDF.h"
#include "endout/endout_util/test/correctIncElasticOutput.h"
#include "endout/endout_util/test/check_MF_Output.h"

using namespace njoy::ENDFtk;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
using Elastic = section::Type<7,2>;

TEST_CASE( "finalizing incoherent elastic scattering data for ENDF" ){
  GIVEN( "alpha Si test case from ENDFB-VIII.0" ){
    //double awr = 27.84423; int za = 147;
    double sigma_b = 2.168771454;
    std::vector<double> temps, debyeWallerVec;
    WHEN( "one temperature" ){
      temps = { 293.6 };
      debyeWallerVec = { 5.14233613 };
      IncoherentElastic mine = writeIncElasticToENDF(sigma_b,temps,debyeWallerVec);
      long lineNumber = 1;
      auto begin = alpha_Si_ENDF_1temp.begin(),
           end   = alpha_Si_ENDF_1temp.end();
      IncoherentElastic good( begin, end, lineNumber, 47, 7, 2 );

      THEN( "the only temperature and only debye waller factor are duplicated" ){
        testIncoherentElasticOutput( mine, good );
      } // THEM
    } // GIVEN
    WHEN( "two temperatures" ){
      temps = { 293.6, 350.0 };
      debyeWallerVec = { 5.14233613, 6.00917141 };

      IncoherentElastic mine = writeIncElasticToENDF(sigma_b,temps,debyeWallerVec);
  
      long lineNumber = 1;
      auto begin = alpha_Si_ENDF_2temps.begin(),
           end   = alpha_Si_ENDF_2temps.end();
      IncoherentElastic good( begin, end, lineNumber, 47, 7, 2 );
      THEN( "temperatures and debye waller factors are written in" ){
        testIncoherentElasticOutput( mine, good );
      } // THEN

    } // WHEN
  } // GIVEN
} // TEST CASE

