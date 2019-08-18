#define CATCH_CONFIG_MAIN
#include <iostream>
#include <vector>
#include "trans/trans.h"
#include "catch.hpp"
#include "generalTools/testing.h"

TEST_CASE( "trans" ){

  std::vector<double> alpha {0.10, 0.20, 0.40, 0.50}, beta {0.15, 0.18, 0.22}, correct;

  double trans_weight = 0.05, delta = 9e-2, c, temp = 296, 
    sc = 1.0, scaling = 1.0, lambda_s = 0.23, tbeta = 0.95, correct_t_eff_val,
    t_eff = 570;
  trans_weight = 0.5;
  tbeta = 0.5;

  std::vector<double> sab {1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,1e-2,1.1e-2,1.2e-2};

  GIVEN( "no diffusion" ){
    c = 0.0;
    /*
    */
    WHEN( "translational weight is small" ){
      trans_weight = 0.05; 
      tbeta = 0.95;
      trans( alpha, beta, trans_weight, delta, c, sc, scaling, lambda_s, tbeta, t_eff, temp,  sab );
      correct = { 1.3599403E+00,  8.4210039E-01,  3.8625375E-01,  1.6486076E+00,  1.3102589E+00, 8.9648052E-01,  1.4734088E+00,  1.3216285E+00,  1.1048539E+00,  1.3628487E+00, 1.2530725E+00,  1.0914529E+00 };
      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
        REQUIRE( 556.3 == Approx(t_eff).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "translational weight is large" ){
      trans_weight = 0.5; 
      tbeta = 0.5;
      trans( alpha, beta, trans_weight, delta, c, sc, scaling, lambda_s, tbeta, t_eff, temp,  sab );
      correct = {1.1704493E+00,  1.1329332E+00,  1.0666042E+00,  8.4594660E-01,  8.3836002E-01, 8.2087238E-01,  5.7370164E-01,  5.7640806E-01,  5.7511395E-01,  4.9844664E-01, 5.0182128E-01,  5.0387717E-01 };
      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
        REQUIRE( 433 == Approx(t_eff).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "alpha and beta values are very small" ){
      alpha = {0.001, 0.002, 0.004, 0.005};
      beta  = {0.0015,0.0018,0.0022};
      trans_weight = 0.05; 
      tbeta = 0.95;
      trans( alpha, beta, trans_weight, delta, c, sc, scaling, lambda_s, tbeta, t_eff, temp,  sab );
      correct = {  3.9432568E+01,  3.9279083E+01,  3.8934953E+01,  2.8027566E+01,  2.7993941E+01,
 2.7848775E+01,  1.9837990E+01,  1.9819911E+01,  1.9795829E+01,  1.7749763E+01,
 1.7735559E+01,  1.7716637E+01 }; 
      THEN( "S(a,b) and effective temperature outputs are correct" ){
        REQUIRE(ranges::equal(sab,correct,equal));
        REQUIRE( 556.3 == Approx(t_eff).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE




/*

TEST_CASE( "trans (old tests)" ){

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
	REQUIRE( correct_t_eff_val == Approx(t_eff).epsilon(1e-6) );
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
*/

