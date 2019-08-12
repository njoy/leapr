#define CATCH_CONFIG_MAIN
#include <iostream>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include "trans/trans.h"
#include "catch.hpp"

auto equal = [](auto x, auto y){return x == Approx(y).epsilon(1e-4);};

TEST_CASE( "trans" ){

  std::vector<double> alpha {0.10, 0.20, 0.40, 0.50}, 
    beta {0.15, 0.18, 0.22}, correct;

  double trans_weight = 0.03, delta = 220.0, diffusion_const = 1.5, temp = 200, 
    sc = 1.0, scaling = 1.0, lambda_s = 0.002, tbeta = 2.1, correct_t_eff_val,
    t_eff = 13.5; ;

  std::vector<double> sab {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2};


  GIVEN( "that the translational motion is diffusive" ){
    WHEN( "temperature is relatively low" ){
      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        lambda_s, tbeta, t_eff, temp,  sab );
      correct = { 0.23049978, 0.25982880, 0.19141505, 0.62197701, 0.58781315, 
        0.39163902, 1.08210491, 0.92354902, 0.64343974, 1.41011128, 1.18123544, 
        0.84745080};
      correct_t_eff_val = 16.12676056;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
	REQUIRE( correct_t_eff_val == Approx(t_eff).epsilon(1e-4) );
      } // THEN
    } // WHEN

    WHEN( "temperature is relatively high" ){

      alpha = {0.8, 1.0, 1.4, 1.5};
      beta = {0.15, 0.19, 0.24, 0.30, 0.31 };
      temp = 800.0;
      t_eff = 117.2;
      std::vector<double> sab { 0.001, 0.002, 0.003, 0.004, 0.006, 0.01, 0.02, 0.03, 0.04, 0.06, 
            0.1, 0.2, 0.3, 0.4, 0.6, 1.1, 1.2, 1.3, 1.4, 1.6 };
      lambda_s = 2.5236078E-3;
      tbeta = 5.1;

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        lambda_s, tbeta, t_eff, temp, sab );
      correct = { 0.92601235894, 0.61808459249, 0.4026621598, 0.2607989292, 
        0.2444620547, 1.0791262271, 0.7478653475, 0.5053651418, 0.3344515293, 
        0.3160592329, 1.3359968514, 1.0435211230, 0.7531824364, 0.5472561336, 
        0.5385673429, 2.0990244422, 1.7440234835, 1.4013032062, 1.0787827232, 
        1.0696952350 };
      correct_t_eff_val = 121.1929824;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
	REQUIRE( correct_t_eff_val == Approx(t_eff).epsilon(1e-4) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "that the translational motion is a free gas" ){
    diffusion_const = 0;

    WHEN( "temperature is relatively low" ){

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        lambda_s, tbeta, t_eff, temp,  sab );
      correct = { 0.92779395, 0.462578, 0.18376845, 1.771662, 1.2519220,
        0.74812670, 2.146880, 1.78983859, 1.384712, 2.2401082, 1.99146590, 
        1.630380 };
      correct_t_eff_val = 16.12676056;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
	REQUIRE( correct_t_eff_val == Approx(t_eff).epsilon(1e-4) );
      } // THEN
    } // WHEN

    WHEN( "temperature is relatively high" ){

      alpha = {0.8, 1.0, 1.4, 1.5};
      beta = {0.15, 0.19, 0.24, 0.30, 0.31 };
      temp = 800.0;
      t_eff = 117.2;
      std::vector<double> sab { 0.001, 0.002, 0.003, 0.004, 0.006, 0.01, 0.02, 
        0.03, 0.04, 0.06, 0.1, 0.2, 0.3, 0.4, 0.6, 1.1, 1.2, 1.3, 1.4, 1.6 };

      lambda_s = 2.5236078E-3;
      tbeta = 5.1;

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        lambda_s, tbeta, t_eff, temp, sab );
      correct = { 1.539329, 1.363224, 1.116673, 0.8210472, 0.7761662, 
        1.446650, 1.318804, 1.1343993, 0.8915002, 0.8546271, 1.362848, 
        1.276393, 1.168248, 0.9998096, 0.9865463, 1.9122958, 1.806748, 
        1.693513, 1.505287, 1.4913750 };
      correct_t_eff_val = 121.19297940773100;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
	REQUIRE( correct_t_eff_val == Approx(t_eff).epsilon(1e-4) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE


