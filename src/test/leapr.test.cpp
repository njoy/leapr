#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "leapr.cpp"

TEST_CASE( "leapr" ){
  int nout, ntempr, iprint, nphon, mat, npr, iel, ncold, nss, nalpha, nbeta, 
      lat, ni, nd, nka, mss;
  double za, awr, spr, aws, sps, delta, trans_weight, c, tbeta, dka, b7;
  std::vector<double> alpha, beta, temp_vec, rho, oscEnergies, oscWeights, 
    kappaVals;
  std::string title;

  nout         = 20;                                                 // Card 1
  title        = "title";                                            // Card 2
  ntempr       = 1;       iprint = 1;     nphon = 3;                 // Card 3
  mat          = 26;      za     = 126.0;                            // Card 4
  awr          = 8.93478; spr    = 6.15;  npr   = 1;   iel = 2;   ncold = 0;
  nss          = 0;       aws    = 0.0;                              // Card 6
  nalpha       = 5;       nbeta  = 5;     lat   = 3;                 // Card 7
  alpha        = { 0.10, 0.20, 0.40, 0.80, 1.60 };                   // Card 8
  beta         = { 0.10, 0.15, 0.30, 0.60, 1.20 };                   // Card 9
  temp_vec     = { 200.0 };                                          // Card 10 
  delta        = 3.8;     ni     = 6;                                // Card 11
  rho          = { 0.002, 0.004, 0.02, 0.04, 0.2, 0.4 };             // Card 12
  trans_weight = 0.3;     c      = 1.0;   tbeta = 20;                // Card 13
  nd           = 2;                                                  // Card 14
  oscEnergies  = { 1.0, 2.0 };                                       // Card 15
  oscWeights   = { 0.3, 0.4 };                                       // Card 16
  nka          = 4;       dka    = 0.01;                             // Card 17
  kappaVals    = { 0.1, 0.2, 0.4, 0.7 };                             // Card 18



  nout         = 24;                                                 // Card 1
  title        = "h in h20, shortened endf model";                   // Card 2
  ntempr       = 1;       iprint = 1;      nphon = 100;              // Card 3
  mat          = 101;     za     = 1001.0;                           // Card 4
  awr          = 0.99917; spr    = 20.449; npr   = 2;   iel = 0;   ncold = 0;
                                                                     // Card 5
  nss          = 1;       b7     = 1.;     aws   = 1.1; sps = 3.8883; mss = 1; 
                                                                     // Card 6
  nalpha       = 5;       nbeta  = 7;      lat   = 1;                // Card 7
  alpha        = { 0.01008, 0.015, 0.0252, 0.033, 0.050406 };        // Card 8
  beta         = { 0.00000, 0.006375, 0.012750, 0.025500, 0.038250, 0.05100, 0.065750 }; // Card 9
  temp_vec     = { 296.0 };                                          // Card 10 
  delta        = 0.00255;     ni     = 67;                           // Card 11
  rho          = { 0.00000, 0.00050, 0.00100, 0.00200, 0.00350, 0.00500, 
                   0.00750, 0.01000, 0.01300, 0.01650, 0.02000, 0.02450, 
                   0.02900, 0.03400, 0.03950, 0.04500, 0.05060, 0.05620, 
                   0.06220, 0.06860, 0.07500, 0.08300, 0.09100, 0.09900, 
                   0.10700, 0.11500, 0.11970, 0.12140, 0.12180, 0.11950, 
                   0.11250, 0.10650, 0.10050, 0.09542, 0.09126, 0.08710, 
                   0.08390, 0.08070, 0.07798, 0.07574, 0.07350, 0.07162, 
                   0.06974, 0.06804, 0.06652, 0.06500, 0.06340, 0.06180, 
                   0.06022, 0.05866, 0.05710, 0.05586, 0.05462, 0.05350, 
                   0.05250, 0.05150, 0.05042, 0.04934, 0.04822, 0.04706, 
                   0.04590, 0.04478, 0.04366, 0.04288, 0.04244, 0.04200, 0. };
                                                                     // Card 12
  trans_weight = 0.055556;     c      = 0.0;    tbeta = 0.444444;    // Card 13
  nd           = 2;                                                  // Card 14
  //oscEnergies  = { 205.0,    0.48};                                  // Card 15
  oscEnergies  = { 35.8,    0.48};                                  // Card 15
  oscWeights   = { 0.166667, 0.333333 };                             // Card 16
  nka          = 4;       dka    = 0.01;                             // Card 17
  kappaVals    = { 0.1, 0.2, 0.4, 0.7 };                             // Card 18

 
  leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, 
      spr, npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, 
      temp_vec, delta, ni, rho, trans_weight, c, tbeta, nd, oscEnergies, 
      oscWeights, nka, dka, kappaVals );





  GIVEN( "inputs" ){
    WHEN( "inputs" ){
      THEN( "results" ){
        REQUIRE( true );
      } // THEN
    } // WHEN
  } // GIVEN 
} // TEST CASE

