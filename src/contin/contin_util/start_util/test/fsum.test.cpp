#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/fsum.h"
#include "range/v3/all.hpp"


TEST_CASE( "fsum" ){
  GIVEN( "a specified beta grid" ){ 
    std::vector<double> p1  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p2 {0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
                        betaGrid (12);
    double tau = 0.5, delta = 2.0;
    for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }

    WHEN( "n = 1 (used for normalizing p(beta)) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        auto beta_P_1 = ranges::view::zip(betaGrid,p1);
        REQUIRE( fsum(1,beta_P_1,tau) == Approx(29610795.32).epsilon(1e-6));

        delta = 1.0;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        auto beta_P_2 = ranges::view::zip(betaGrid,p2);
        REQUIRE( fsum(1,beta_P_2,tau) == Approx(39.387006).epsilon(1e-6) );

      } // THEN
    } // WHEN
    for ( auto& entry : p2 ){ entry *= 0.1; } 
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        auto beta_P_1 = ranges::view::zip(betaGrid,p1);
        REQUIRE( fsum(0,beta_P_1,tau) == Approx(1444532.8400).epsilon(1e-6) );

        delta = 0.1;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        auto beta_P_2 = ranges::view::zip(betaGrid,p2);
        REQUIRE( fsum(0,beta_P_2,tau) == Approx(0.035514341).epsilon(1e-6) );
  
      } // THEN
    } // WHEN
    WHEN( "n = 2 (used for effective temperature calculation) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        auto beta_P_1 = ranges::view::zip(betaGrid,p1);
        REQUIRE( fsum(2,beta_P_1,tau) == Approx(612298146.17046).epsilon(1e-6) );
 
        delta = 0.7;
        for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
        auto beta_P_2 = ranges::view::zip(betaGrid,p2);
        REQUIRE( fsum(2,beta_P_2,tau) == Approx(3.21945524).epsilon(1e-6) );
  
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
