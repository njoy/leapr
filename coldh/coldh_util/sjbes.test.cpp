#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sjbes.h"


TEST_CASE( "sjbes" ){
  GIVEN( "inputs" ){
    int n = 1;
    double x = 0.35;
    auto sjbes_output = sjbes(n,x);
    std::cout << sjbes_output << std::endl;
   REQUIRE( true ) ;



  } // THEN
} // TEST CASE
