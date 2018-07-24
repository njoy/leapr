#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/fsum.h"


TEST_CASE( "fsum" ){
  GIVEN( "input value p in the form of a ranges view" ){ 

    double tau = 0.5, delta = 2.0;
    auto p = ranges::view::iota(1,13);

    WHEN( "n = 1 (used for normalizing p(beta)) " ){

      auto p1 = p | ranges::view::take_exactly(6) 
                  | ranges::view::transform([](auto x){ return 0.1*x; });

      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum(1,p,tau,delta) == Approx(29610795.32).epsilon(1e-6));
        REQUIRE( fsum(1,p1,tau,1.0) == Approx(39.387006).epsilon(1e-6) );

      } // THEN
    } // WHEN

    auto p1 = p | ranges::view::take_exactly(6) 
                | ranges::view::transform([](auto x){ return 0.01*x; });


    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        REQUIRE( fsum(0,p,tau,delta) == Approx(1444532.8400).epsilon(1e-6) );
        REQUIRE( fsum(0,p1,tau,0.1) == Approx(0.035514341).epsilon(1e-6) );

 
      } // THEN
    } // WHEN
    WHEN( "n = 2 (used for effective temperature calculation) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum(2,p,tau,delta) == Approx(612298146.17046).epsilon(1e-6) );
        REQUIRE( fsum(2,p1,tau,0.7) == Approx(3.21945524).epsilon(1e-6) );
  
 
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "input value p in the form of a vector" ){ 
    std::vector<double> p  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p1 {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    double tau = 0.5, delta = 2.0;

    WHEN( "n = 1 (used for normalizing p(beta)) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum(1,p,tau,delta) == Approx(29610795.32).epsilon(1e-6));
        REQUIRE( fsum(1,p1,tau,1.0) == Approx(39.387006).epsilon(1e-6) );

      } // THEN
    } // WHEN
    for ( auto& entry : p1 ){ entry *= 0.1; } 
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum(0,p,tau,delta) == Approx(1444532.8400).epsilon(1e-6) );
        REQUIRE( fsum(0,p1,tau,0.1) == Approx(0.035514341).epsilon(1e-6) );

  
      } // THEN
    } // WHEN
    WHEN( "n = 2 (used for effective temperature calculation) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){

        REQUIRE( fsum(2,p,tau,delta) == Approx(612298146.17046).epsilon(1e-6) );
        REQUIRE( fsum(2,p1,tau,0.7) == Approx(3.21945524).epsilon(1e-6) );

  
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
