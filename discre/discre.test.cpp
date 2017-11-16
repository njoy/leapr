#define CATCH_CONFIG_MAIN
#include "../catch.hpp" 
#include "discre.h"
#include <iostream> 

TEST_CASE( "discre" ){
  GIVEN( "input" ){
    std::vector<double> alpha { 0.1, 0.2, 0.4, 0.8, 1.6 };
    std::vector<double> beta  { 0.10, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> osc_energies { 0.035, 0.05 };
    std::vector<double> osc_weights  { 0.2, 0.8 };
    std::vector<double> t_eff_vec { 81178.935219 };
    std::vector<double> temp_vec { 200.0 };
    std::vector<std::vector<std::vector<double>>> sym_sab( alpha.size(),
      std::vector<std::vector<double>> ( beta.size(), 
        std::vector<double> (1, 0.0) ) );
    int i = 1;
    for ( auto a = 0; a < alpha.size(); ++a ){
      for ( auto b = 0; b < beta.size(); ++b ){
        sym_sab[a][b][0] = i;
        i += 1;
      }
    }

    double lambda_s = 2.2941534E-3;
    double sc = 1.0;
    double scaling = 1.0;
    double tev = 1.723477E-2;
    double tbeta = 2.0;
    double twt = 0.3;
    int itemp = 0;
    discre( sc, scaling, alpha, beta, tev, lambda_s, osc_energies, osc_weights, 
        tbeta, t_eff_vec, temp_vec, itemp, sym_sab, twt );

  } // GIVEN
} // TEST CASE
