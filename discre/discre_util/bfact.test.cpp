#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "bfact.h"

TEST_CASE( "bfact" ){
  GIVEN( "inputs" ){
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    double x   = 8.213274e-3;
    double dwc = 1.282379e-2;
    double beta_i = 2.030778478;
    bfact(x, dwc, beta_i );
  } // GIVEN
} // TEST CASE
