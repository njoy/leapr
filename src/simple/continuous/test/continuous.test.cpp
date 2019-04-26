#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/continuous/continuous.h"




TEST_CASE( "continuous" ){
  GIVEN( "vectors" ){

    std::vector<double> rhoEnergy { 0.0, 1.0, 2.0 }, 
                        rho       { 1.2, 1.4, 1.0 },
                        alphas    { 0.01, 0.02, 0.03, 0.04},
                        betas     { 0.0, 0.5, 1.0};
    //continuous(rhoEnergy,rho,alphas,betas);
    REQUIRE(true);
  } // GIVEN
  /*

  GIVEN( "int" ){

    std::vector<double> rhoEnergy { 0.00, 0.01, 0.02, 0.04, 0.08 }, 
                        rho       { 1.00, 2.00, 1.50, 3.00, 0.00 },
                        alphas    { 0.01, 0.02, 0.03, 0.04 },
                        betas     { 0.0, 0.1, 0.7, 1.0, 2.4 };
    double kbT = 0.025;
    continuous(rhoEnergy,rho,alphas,betas,kbT);
    REQUIRE(true);
  } // GIVEN
  */

} // TEST



