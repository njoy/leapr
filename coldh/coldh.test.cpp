#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "coldh.h"

TEST_CASE( "coldh" ){
  GIVEN( "inputs" ){
    int itemp = 0;
    double temp = 200.0;
    double tev = 1.723477e-2;
    double sc = 1.0;
    int ncold = 2;
    coldh( itemp, temp, tev, sc, ncold );


  } // GIVEN
} // TEST CASE
