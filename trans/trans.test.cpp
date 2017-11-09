#define CATCH_CONFIG_MAIN
#include <iostream>
#include <vector>
#include "trans.h"
#include "../catch.hpp"

void equal( double a, double b ){
    if( b == 0 ){ REQUIRE( (a-b) < 1e-4 ); }
    if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-4 ); }
}

TEST_CASE( "trans" ){
  GIVEN( "inputs" ){
    std::vector<double> alpha {0.10,  0.20,  0.40, 0.50};
    std::vector<double> beta  {0.15,  0.18,  0.22};
    std::vector<double> temps {200.0}; 
    std::vector<std::vector<std::vector<double>>> sym_sab ( alpha.size(),
      std::vector<std::vector<double>> ( beta.size(),
        std::vector<double> ( temps.size(), 0.0 ) ) );
    sym_sab = { { { 0.1 }, { 0.2 }, { 0.3 } },
                { { 0.4 }, { 0.5 }, { 0.6 } },
                { { 0.7 }, { 0.8 }, { 0.9 } },
                { { 1.0 }, { 1.1 }, { 1.2 } } };
    double trans_weight = 0.03;
    double delta = 220.0;
    double diffusion_const = 1.5;
    double sc = 1.0;
    double scaling = 1.0;
    int itemp = 0;
    double lambda_s = 0.002;
    std::vector<double> t_eff_vec = {13.5};
    std::vector<double> temp_vec = {200.0};
    double tbeta = 2.1;

    trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, itemp, 
      lambda_s, tbeta, t_eff_vec, temp_vec,  sym_sab );
    std::vector<double> correct{ 0.23049978, 0.25982880, 0.19141505, 
      0.62197701, 0.58781315, 0.39163902, 1.08210491, 0.92354902, 
      0.64343974, 1.41011128, 1.18123544, 0.84745080};
    double correct_t_eff_val = 16.12676056;
    int i = 0;
    for ( auto beta_temp : sym_sab ){
      for ( auto temp_vec : beta_temp ){
        equal( temp_vec[0], correct[i] );
        i += 1;
      } // temp
    } // beta and temp
    equal( t_eff_vec[0], correct_t_eff_val );
  } // GIVEN
  GIVEN( "other input" ){
    std::vector<double> alpha {0.8, 1.0, 1.4, 1.5};
    std::vector<double> beta {0.15, 0.19, 0.24, 0.30, 0.31 };
    std::vector<double> temps {600.0}; 
    std::vector<std::vector<std::vector<double>>> sym_sab ( alpha.size(),
    std::vector<std::vector<double>> ( beta.size(),
       std::vector<double> ( temps.size(), 0.0 ) ) );
    sym_sab = { { { 0.001 }, { 0.002 }, { 0.003 }, { 0.004 }, { 0.006 } },
                { { 0.01  }, { 0.02  }, { 0.03  }, { 0.04  }, { 0.06  } },
                { { 0.1   }, { 0.2   }, { 0.3   }, { 0.4   }, { 0.6   } },
                { { 1.1   }, { 1.2   }, { 1.3   }, { 1.4   }, { 1.6 } } };
    double trans_weight = 0.03;
    double delta = 220.0;
    double diffusion_const = 1.5;
    double sc = 1.0;
    double scaling = 1.0;
    int itemp = 0;
    double lambda_s = 2.5236078E-3;
    std::vector<double> t_eff_vec = {117.2};
    std::vector<double> temp_vec = {800.0};
    double tbeta = 5.1;

    trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, itemp, 
      lambda_s, tbeta, t_eff_vec, temp_vec, sym_sab );
    std::vector<double> correct {0.92601235894, 0.61808459249, 
      0.4026621598, 0.2607989292, 0.2444620547, 1.0791262271, 
      0.7478653475, 0.5053651418, 0.3344515293, 0.3160592329, 
      1.3359968514, 1.0435211230, 0.7531824364, 0.5472561336, 
      0.5385673429, 2.0990244422, 1.7440234835, 1.4013032062, 
      1.0787827232, 1.0696952350 };
    double correct_t_eff_val = 121.1929824;
    int i = 0;
    for ( auto beta_temp_vec : sym_sab ){
      for ( auto temp_vec : beta_temp_vec ){
        equal( temp_vec[0], correct[i] );
        i += 1;
      } // for temp value
    } // for beta, temp values
    equal( t_eff_vec[0], correct_t_eff_val );
  } // GIVEN
} // TEST CASE
 
