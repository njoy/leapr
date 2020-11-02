#include "catch.hpp"
#include "continuous/continuous_util/start.h"
#include "generalTools/testing.h"


TEST_CASE( "getting Debye Waller coefficient" ){
  GIVEN( "a specified beta grid" ){ 
    std::vector<double> p1  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p2  {0.01, 0.02, 0.03, 0.04, 0.05, 0.06},
                        betaGrid1 = makeGrid(int(p1.size()),2.0),
                        betaGrid2 = makeGrid(int(p2.size()),0.1);
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        auto beta_P_1 = ranges::view::zip(betaGrid1,p1);
        auto beta_P_2 = ranges::view::zip(betaGrid2,p2);
        REQUIRE( getDebyeWaller(beta_P_1) == Approx(1444532.840).epsilon(1e-6) );
        REQUIRE( getDebyeWaller(beta_P_2) == Approx(0.035514341).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE





TEST_CASE( "getting effective temperature" ){
  GIVEN( "a specified beta grid" ){ 
    std::vector<double> p1  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p2  {0.01, 0.02, 0.03, 0.04, 0.05, 0.06},
                        betaGrid1 = makeGrid(int(p1.size()),2.0),
                        betaGrid2 = makeGrid(int(p2.size()),0.7);
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        auto beta_P_1 = ranges::view::zip(betaGrid1,p1);
        auto beta_P_2 = ranges::view::zip(betaGrid2,p2);
        REQUIRE( getEffectiveTemp(beta_P_1) == Approx(612298146.17046).epsilon(1e-6) );
        REQUIRE( getEffectiveTemp(beta_P_2) == Approx(3.21945524).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE







TEST_CASE( "normalize" ){
  std::vector<double> p, correct;
  double continWgt, delta;
  GIVEN( "some vector p, spacing delta, normalizing value continWgt" ){
    p = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
    continWgt = 0.8, delta = 0.5;
    std::vector<double> betaGrid = makeGrid(int(p.size()),delta);

    auto beta_P_zipped = ranges::view::zip(betaGrid,p);
    auto normalizedP = normalize( beta_P_zipped, continWgt );
    correct= {2.62155807E-2, 5.24311614E-2, 7.86467421E-2, 0.104862322, 
              0.1310779, 0.15729348};

    THEN( "each value in vector has as most a 1e-6 percent error" ){
      REQUIRE(ranges::equal(normalizedP, correct, equal));
    } // THEN

    {
      p = {7.802837847e-3, 1.560567569e-2, 2.340851354e-2, 3.121135138e-2, 
        3.901418923e-2, 4.681702708e-2, 5.461986492e-2, 6.242270277e-2, 
        7.022554062e-2, 7.802837847e-2, 8.583121631e-2, 9.363405416e-2};

      delta = 0.2; continWgt = 2.0;
      betaGrid.resize(p.size());
      for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
    
      auto beta_P_zipped = ranges::view::zip(betaGrid,p);
      auto normalizedP = normalize( beta_P_zipped, continWgt );
      correct = {5.2978935E-2, 0.10595787, 0.15893680, 0.21191574, 
        0.26489467, 0.31787361, 0.37085254, 0.42383148, 
        0.47681042, 0.52978935, 0.58276829, 0.63574722};


      THEN( "each value in vector has as most a 1e-6 percent error" ){
        REQUIRE(ranges::equal(normalizedP, correct, equal));
      } // THEN
    }
  } // GIVEN


  GIVEN( "two vectors where one is the first scaled by some constant" ){
    std::vector<double> p1 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12}, 
                        p2 = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    continWgt = 5.0, delta = 0.8;
    auto betaGrid = ranges::view::iota(0,int(p1.size())) 
             | ranges::view::transform([delta](auto x){ return x * delta;});

    auto normalizedP1 = normalize( ranges::view::zip(betaGrid,p1), continWgt );
    auto normalizedP2 = normalize( ranges::view::zip(betaGrid,p2), continWgt );

    std::vector<double> correct= { 3.095826174E-2, 6.191652348E-2, 
      9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};

    THEN( "outputs should be same, and both within 1e-6 to true value" ){
      REQUIRE(ranges::equal(normalizedP1, correct, equal));
      REQUIRE(ranges::equal(normalizedP2, correct, equal));
    } // THEN
  } // GIVEN
} // TEST CASE











TEST_CASE( "start function" ){

  std::vector<double> p, correct;
  double delta, tev, tbeta, correct_lambda_s, correct_t_eff;

  GIVEN( "phonon distribution, spacing, and temp in ev" ){

    p = {0.1, 0.2, 0.3, 0.5, 0.8, 1.3};
    delta = 0.0001; tev = 0.001; tbeta = 1.0;
    std::vector<double> betaGrid = makeGrid(int(p.size()), delta/tev);
    auto output = start(p, tbeta, betaGrid);
    correct = {1.9662112, 2.06616, 0.8135182, 0.632185, 0.5964, 0.649624};
    correct_lambda_s = 41.517752; correct_t_eff = 1.0118507;

    THEN( "T1, debye waller, and effective temp are correctly computed" ){
      REQUIRE( std::get<0>(output) == Approx(correct_lambda_s).epsilon(1e-6) );
      REQUIRE( std::get<1>(output) == Approx(correct_t_eff).epsilon(1e-6) );
      REQUIRE( ranges::equal(correct,std::get<2>(output),equal) );
    } // THEN

  } // GIVEN

  GIVEN( "more frequency distribution values" ){

    p = {.001, .003, .006, .008, .01, .03, .06, .08, 1.0, 3.0, 6.0, 8.0};
    tbeta = 2.0;

    WHEN( "spacing and temperature are very medium" ){
      delta = 0.015, tev = 1.723477E-2;
      std::vector<double> betaGrid = makeGrid(int(p.size()), delta/tev);

      auto output = start(p, tbeta, betaGrid);

      correct = {2.691905E-3, 4.031148E-3, 2.841214E-3, 2.247653E-3, 
                 2.014354E-3, 4.746877E-3, 7.851893E-3, 8.945384E-3, 
                 9.771149E-2, 0.2604205, 0.4686490, 0.5680046};
      correct_lambda_s = 0.238142827, correct_t_eff = 4.298571365;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        REQUIRE( std::get<0>(output) == Approx(correct_lambda_s).epsilon(1e-6) );
        REQUIRE( std::get<1>(output) == Approx(correct_t_eff).epsilon(1e-6) );
        REQUIRE( ranges::equal(correct,std::get<2>(output),equal) );
      } // THEN
    } // WHEN

    WHEN( "temperature is changed" ){
      delta = 0.015; tev = 4.3086925E-2;
      std::vector<double> betaGrid = makeGrid(int(p.size()), delta/tev);
      auto output = start(p, tbeta, betaGrid);
      correct = {1.54443108E-2, 1.82883208E-2, 1.07199712E-2, 7.37431363E-3, 
                 5.96172958E-3, 1.30407541E-2, 2.04553735E-2, 2.24449206E-2, 
                 0.238766209, 0.62462916, 1.10947058, 1.332378};
      correct_lambda_s = 0.6485580; correct_t_eff = 1.836531025;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        REQUIRE( std::get<0>(output) == Approx(correct_lambda_s).epsilon(1e-6) );
        REQUIRE( std::get<1>(output) == Approx(correct_t_eff).epsilon(1e-6) );
        REQUIRE( ranges::equal(correct,std::get<2>(output),equal) );
      } // THEN
    } // WHEN

    WHEN( "spacing is very large" ){
      tev = 4.3086925E-2; delta = 4.015;
      std::vector<double> betaGrid = makeGrid(int(p.size()), delta/tev);
      auto output = start(p, tbeta, betaGrid);

      correct = {2.36543334e-7, 2.20419880e-5, 2.20419880e-5, 1.95928782e-5, 
                 1.83683233e-5, 4.40839761e-5, 7.34732935e-5, 8.39694783e-5, 
                 9.18416169e-4, 2.44910978e-3, 4.40839761e-3, 5.34351226e-3};
      correct_lambda_s = 2.2081257E-3; correct_t_eff = 459.94244303;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        REQUIRE( std::get<0>(output) == Approx(correct_lambda_s).epsilon(1e-6) );
        REQUIRE( std::get<1>(output) == Approx(correct_t_eff).epsilon(1e-6) );
        REQUIRE( ranges::equal(correct,std::get<2>(output),equal) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST_CASE
