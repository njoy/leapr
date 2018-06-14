#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coh/coh.h"


TEST_CASE( "coh" ){
  GIVEN( "" ){
    WHEN( "" ){
      double temp = 296.0, emax = 1.2;
      int lat = 1, ne = 145, nex = 2, natom = 1;
      std::vector<double> fl, bufo(1000,0.0), bufn(1000,0.0);
        bufn[0] = 1e-5;      bufn[1] = 78.42469; bufn[2] = 1.0625e-5; 
        bufn[3] = 76.09049;  bufn[4] = 1.125e-5; bufn[5] = 73.95383; 
        bufn[6] = 1.1875e-5; bufn[7] = 71.98835; bufn[8] = 1.25e-5;
      std::fstream iold("test_coh_1a"), inew("test_coh_1b");
      coh( lat, inew, ne, nex, temp, iold, emax, natom, fl, bufo, bufn );
      REQUIRE( true );
      
    } // WHEN
  } // GIVEN
} // TEST CASE
