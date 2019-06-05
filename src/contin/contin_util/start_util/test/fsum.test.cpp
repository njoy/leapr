#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/fsum.h"


TEST_CASE( "fsum3" ){
  GIVEN( "a specified beta grid" ){ 
    std::vector<double> p  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p1 {0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
                        betaGrid (12);
    double tau = 0.5, delta = 2.0;
    for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }

    WHEN( "n = 1 (used for normalizing p(beta)) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        REQUIRE( fsum3(1,p,tau,betaGrid) == Approx(29610795.32).epsilon(1e-6));

        delta = 1.0;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        REQUIRE( fsum3(1,p1,tau,betaGrid) == Approx(39.387006).epsilon(1e-6) );

      } // THEN
    } // WHEN
    for ( auto& entry : p1 ){ entry *= 0.1; } 
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum3(0,p,tau,betaGrid) == Approx(1444532.8400).epsilon(1e-6) );

        delta = 0.1;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        REQUIRE( fsum3(0,p1,tau,betaGrid) == Approx(0.035514341).epsilon(1e-6) );
  
      } // THEN
    } // WHEN
    WHEN( "n = 2 (used for effective temperature calculation) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum3(2,p,tau,betaGrid) == Approx(612298146.17046).epsilon(1e-6) );
 
        delta = 0.7;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        REQUIRE( fsum3(2,p1,tau,betaGrid) == Approx(3.21945524).epsilon(1e-6) );
  
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
