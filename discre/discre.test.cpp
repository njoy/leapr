#define CATCH_CONFIG_MAIN
#include "../catch.hpp" 
#include "discre.h"
#include <iostream> 

TEST_CASE( "discre" ){
  GIVEN( "input" ){
    std::vector<double> alpha { 0.1, 0.2, 0.4, 0.8, 1.6 };
    std::vector<double> beta  { 0.0, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> osc_energies { 0.035, 0.05 };
    std::vector<double> osc_weights  { 0.2, 0.8 };
    double lambda_s = 2.2941534E-3;
    double sc = 1.0;
    double tev = 1.723477E-2;

    discre( sc, alpha, beta, tev, lambda_s, osc_energies, osc_weights );

    std::cout << "Hello, world" << std::endl;
  } // GIVEN
} // TEST CASE
