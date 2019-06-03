#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/continuous/continuous.h"


TEST_CASE( "continuous" ){
  GIVEN( "vectors" ){

    std::vector<double> rhoEnergy { 0.0, 0.01, 0.02, 0.03 }, 
                        rho       { 0.0, 1.2, 1.4, 1.0 },
                        alphas    { 0.1, 0.2, 0.3, 0.4},
                        betas     { 0.0, 0.5, 1.0, 5.0};
    betas = rho;
    int lat = 1;
    double kbT = 296*8.617333262145E-5, tbeta = 1.0, arat = 1.0;
    continuous(rhoEnergy,rho,alphas,betas,kbT,tbeta,lat,arat);
    REQUIRE(true);
  } // GIVEN
} // TEST



