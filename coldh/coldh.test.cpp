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
    double tbeta = 2.0;
    double trans_weight = 0.3;
    double scaling = 1.0;
    std::vector<double> tempf {193093.99765};
    std::vector<double> alpha {0.1, 0.2, 0.4, 0.8, 1.6};
    std::vector<double> tempr {200.0};
    coldh( itemp, temp, tev, sc, ncold, trans_weight, tbeta, tempf, tempr, 
      scaling, alpha );


  } // GIVEN
} // TEST CASE
