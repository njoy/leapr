#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/continuous/continuous.h"


TEST_CASE( "continuous" ){
  GIVEN( "int" ){

    std::vector<double> rhoEnergy { 0.0, 1.0, 2.0 }, 
                        rho       { 1.2, 1.4, 1.0 },
                        alphas    { 0.01, 0.02, 0.03, 0.04},
                        betas     { 0.0, 0.5, 1.0};
    continuous(rhoEnergy,rho,alphas,betas);
    REQUIRE(true);
  } // GIVEN
} // TEST



