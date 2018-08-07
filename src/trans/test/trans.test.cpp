#define CATCH_CONFIG_MAIN
#include <iostream>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include "trans/trans.h"
#include "catch.hpp"


void equalSAB( const Eigen::Tensor<double,3>& sab,
  const std::vector<double>& correct ){

  REQUIRE( sab.dimension(0)*sab.dimension(1)*sab.dimension(2) == correct.size() );

  int l = 0;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      for ( int k = 0; k < sab.dimension(2); ++k ){
        REQUIRE( sab(i,j,k) == Approx(correct[l]).epsilon(1e-4) );
	l += 1;
      }
    }
  }
}


TEST_CASE( "trans" ){

  std::vector<double> alpha {0.10, 0.20, 0.40, 0.50}, 
    beta {0.15, 0.18, 0.22}, temps {200.0}, t_eff_vec = {13.5}, 
    correct;

  double trans_weight = 0.03, delta = 220.0, diffusion_const = 1.5, 
    sc = 1.0, scaling = 1.0, lambda_s = 0.002, tbeta = 2.1, correct_t_eff_val;

  Eigen::Tensor<double,3> sym_sab( int(alpha.size()), int(beta.size()), int(temps.size()) );
  sym_sab.setValues ( { { { 0.1 }, { 0.2 }, { 0.3 } },
  		        { { 0.4 }, { 0.5 }, { 0.6 } },
                        { { 0.7 }, { 0.8 }, { 0.9 } },
                        { { 1.0 }, { 1.1 }, { 1.2 } } } );

  int itemp = 0;

  /*

  GIVEN( "that the translational motion is diffusive" ){
    WHEN( "temperature is relatively low" ){

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        itemp, lambda_s, tbeta, t_eff_vec, temps,  sym_sab );
      correct = { 0.23049978, 0.25982880, 0.19141505, 0.62197701, 0.58781315, 
        0.39163902, 1.08210491, 0.92354902, 0.64343974, 1.41011128, 1.18123544, 
        0.84745080};
      correct_t_eff_val = 16.12676056;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        equalSAB( sym_sab, correct );
	REQUIRE( correct_t_eff_val == Approx(t_eff_vec[0]).epsilon(1e-4) );
      } // THEN
    } // WHEN

    WHEN( "temperature is relatively high" ){

      alpha = {0.8, 1.0, 1.4, 1.5};
      beta = {0.15, 0.19, 0.24, 0.30, 0.31 };
      Eigen::Tensor<double,3> sym_sab( int(alpha.size()), int(beta.size()), int(temps.size()) );
      temps = {800.0};
      t_eff_vec = {117.2};
      sym_sab.setValues( { { {0.001}, {0.002}, {0.003}, {0.004}, {0.006} },
                           { {0.01 }, {0.02 }, {0.03 }, {0.04 }, {0.06 } },
                           { {0.1  }, {0.2  }, {0.3  }, {0.4  }, {0.6  } },
                           { {1.1  }, {1.2  }, {1.3  }, {1.4  }, {1.6  } } } );
      lambda_s = 2.5236078E-3;
      tbeta = 5.1;

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        itemp, lambda_s, tbeta, t_eff_vec, temps, sym_sab );
      correct = { 0.92601235894, 0.61808459249, 0.4026621598, 0.2607989292, 
        0.2444620547, 1.0791262271, 0.7478653475, 0.5053651418, 0.3344515293, 
        0.3160592329, 1.3359968514, 1.0435211230, 0.7531824364, 0.5472561336, 
        0.5385673429, 2.0990244422, 1.7440234835, 1.4013032062, 1.0787827232, 
        1.0696952350 };
      correct_t_eff_val = 121.1929824;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        equalSAB( sym_sab, correct );
	REQUIRE( correct_t_eff_val == Approx(t_eff_vec[0]).epsilon(1e-4) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "that the translational motion is a free gas" ){
    diffusion_const = 0;

    WHEN( "temperature is relatively low" ){

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        itemp, lambda_s, tbeta, t_eff_vec, temps,  sym_sab );
      correct = { 0.92779395, 0.462578, 0.18376845, 1.771662, 1.2519220,
        0.74812670, 2.146880, 1.78983859, 1.384712, 2.2401082, 1.99146590, 
        1.630380 };
      correct_t_eff_val = 16.12676056;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        equalSAB( sym_sab, correct );
	REQUIRE( correct_t_eff_val == Approx(t_eff_vec[0]).epsilon(1e-4) );
      } // THEN
    } // WHEN

    WHEN( "temperature is relatively high" ){

      alpha = {0.8, 1.0, 1.4, 1.5};
      beta = {0.15, 0.19, 0.24, 0.30, 0.31 };
      temps = {800.0};
      t_eff_vec = {117.2};
      Eigen::Tensor<double,3> sym_sab( int(alpha.size()), int(beta.size()), int(temps.size()) );
      sym_sab.setValues( { { {0.001}, {0.002}, {0.003}, {0.004}, {0.006} },
                           { {0.01 }, {0.02 }, {0.03 }, {0.04 }, {0.06 } },
                           { {0.1  }, {0.2  }, {0.3  }, {0.4  }, {0.6  } },
                           { {1.1  }, {1.2  }, {1.3  }, {1.4  }, {1.6  } } } );
      lambda_s = 2.5236078E-3;
      tbeta = 5.1;

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        itemp, lambda_s, tbeta, t_eff_vec, temps, sym_sab );
      correct = { 1.539329, 1.363224, 1.116673, 0.8210472, 0.7761662, 
        1.446650, 1.318804, 1.1343993, 0.8915002, 0.8546271, 1.362848, 
        1.276393, 1.168248, 0.9998096, 0.9865463, 1.9122958, 1.806748, 
        1.693513, 1.505287, 1.4913750 };
      correct_t_eff_val = 121.19297940773100;

      THEN( "S(a,b) and effective temperature outputs are correct" ){
        equalSAB( sym_sab, correct );
	REQUIRE( correct_t_eff_val == Approx(t_eff_vec[0]).epsilon(1e-4) );
      } // THEN
    } // WHEN
  } // GIVEN
  */
  GIVEN( "H in H2O inputs" ){
    diffusion_const = 0;
    std::vector<double> alpha { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 0.100812, 
    0.151218, 0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 
    0.504060, 0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 
    0.976435, 1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 
    1.873790, 2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 
    3.578320, 3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 
    6.816500, 7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 
    12.96740, 14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 
    24.65260, 27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 
    46.85030, 50.0 },
  beta { 0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, 
    0.0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 
    0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 
    1.129090, 1.209740, 1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 
    1.693640, 1.774290, 1.854940, 1.935590, 2.016240, 2.096890, 2.177540, 
    2.258190, 2.338840, 2.419490, 2.500140, 2.580790, 2.669500, 2.767090, 
    2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 3.595360, 3.785490, 
    3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 5.769970, 
    6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 
    10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 17.17330, 
    18.72180, 20.42450, 22.29760, 24.35720, 25.0 };

  Eigen::Tensor<double,3> sym_sab( int(alpha.size()), int(beta.size()), int(temps.size()) );
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sym_sab(i,j,0) = (i+j+2)*0.0005;
    } 
  }
  std::cout << beta.size() << std::endl;
  std::cout << sym_sab(49,65,0) << std::endl;
   scaling = 0.99186670867058835; 
   sc = 0.99186670867058835;

    trans_weight = 5.5556e-2; diffusion_const = 0.0;
    lambda_s = 0.23520650571218535; tbeta = 0.444444;

    temps = {296};

    WHEN( "temperature is relatively low" ){

      trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling, 
        itemp, lambda_s, tbeta, t_eff_vec, temps,  sym_sab );
      std::vector<double> correct_alpha_0_49_beta_0_49 = { 11.942047, 9.6832173, 
        7.3444064, 6.0882386, 4.7270718, 3.7768044, 3.1493396, 2.5589613, 
        1.9599579, 1.5438730, 0.97229950, 0.59355179, 0.35779156, 0.21598531, 
        0.13219498, 8.3561190e-2, 5.6782269e-2, 4.2119690e-2, 3.4104938e-2, 
        2.9790332e-2, 2.7576664e-2, 2.6581416e-2, 2.6286551e-2, 2.6432962e-2, 
        2.6842441e-2, 2.7429248e-2, 2.8129284e-2, 2.8902840e-2, 2.9728511e-2, 
        3.0590201e-2, 3.1477126e-2, 3.2380600e-2, 3.3291412e-2, 3.4207080e-2, 
        3.5122241e-2, 3.6034452e-2, 3.6938781e-2, 3.7826413e-2, 3.8699194e-2, 
        3.9551834e-2, 4.0361894e-2, 4.1135875e-2, 4.1877405e-2, 4.2595236e-2, 
        4.3300122e-2, 4.4007130e-2, 4.4730310e-2, 4.5481599e-2, 4.6271257e-2, 
        4.7105334e-2 };
      //correct_t_eff_val = 16.12676056;
      //std::cout << std::setprecision(15) << sym_sab(0,0,0) << std::endl;

      //std::cout << sym_sab.dimension(0) << std::endl;
      //std::cout << sym_sab(48,48,0) << std::endl;
      //std::cout << sym_sab(49,49,0) << std::endl;
      THEN( "S(a,b) and effective temperature outputs are correct" ){
        for ( size_t i = 0; i < correct_alpha_0_49_beta_0_49.size(); ++i ){ 
          REQUIRE( correct_alpha_0_49_beta_0_49[i] == Approx(sym_sab(i,i,0)).epsilon(1e-3) );
        }
        //equalSAB( sym_sab, correct );
	//REQUIRE( correct_t_eff_val == Approx(t_eff_vec[0]).epsilon(1e-4) );
      } // THEN
    } // WHEN
  } // GIVEN

} // TEST CASE
 
