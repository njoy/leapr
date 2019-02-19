#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/simpleLEAPR.h"


TEST_CASE( "leapr" ){
  GIVEN( "int" ){
    std::vector<double> rho = { 1.0,2.0,3.0 }, alphas {0.1,0.2}, betas {0.0, 1.0};
    simpleLEAPR(rho,alphas,betas);

    REQUIRE(true);
  } // GIVEN
} // TEST



