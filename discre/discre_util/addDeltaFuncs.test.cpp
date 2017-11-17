#define CATCH_CONFIG_MAIN
#include "../../catch.hpp" 
#include "addDeltaFuncs.h"

TEST_CASE( "helper function to add delta functions to scattering law" ){
  GIVEN( "inputs" ){
    double twt = -0.1;
    double dwf = 0.9;
    int n = 1;
    std::vector<double> bes ( 500, 0.0 ); 
    std::vector<double> betan { 0.1, 0.2, 0.3, 0.5 };
    std::vector<double> wts { 0.1, 0.2, 0.3, 0.5 };
    std::vector<double> sexpb { 0.1, 0.2, 0.3, 0.5 };
    double val = 9; int sign = -1;
    for ( auto i = 1; i < 20; ++i ){ 
      bes[i] = sign * val;
      sign *= -1; 
      val -= 1;
    }
    for ( auto entry : bes ){ std::cout << entry << std::endl; }
    addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
  } // GIVEN
} // TEST CASE
