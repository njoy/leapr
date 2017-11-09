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
    std::vector<double> t_eff_vec { 81178.935219 };
    std::vector<double> temp_vec { 200.0 };
    double lambda_s = 2.2941534E-3;
    double sc = 1.0;
    double scaling = 1.0;
    double tev = 1.723477E-2;
    double tbeta = 2.0;
    int itemp = 0;
    discre( sc, scaling, alpha, beta, tev, lambda_s, osc_energies, osc_weights, 
        tbeta, t_eff_vec, temp_vec, itemp );

  } // GIVEN
} // TEST CASE
