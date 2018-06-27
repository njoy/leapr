#include "catch.hpp"
#include "contin/contin_util/start.h"
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>

void check( const std::vector<double>& p, const std::vector<double>& correct,
  const std::tuple<double,double>& output, const double& lambda, 
  const double& effectiveTemp ){

  REQUIRE( p.size() == correct.size() );
  for ( size_t i = 0; i < p.size(); ++i ){
    REQUIRE( p[i] == Approx(correct[i]).epsilon(1e-6) );
  }
  REQUIRE( std::get<0>(output) == Approx(lambda).epsilon(1e-6) );
  REQUIRE( std::get<1>(output) == Approx(effectiveTemp).epsilon(1e-6) );

}


TEST_CASE( "start function" ){

  std::vector<double> p, correct;
  std::tuple<double,double> output;
  double delta, tev, tbeta, correct_lambda_s, correct_t_eff;

  GIVEN( "phonon distribution, spacing, and temp in ev" ){

    p = {0.1, 0.2, 0.3, 0.5, 0.8, 1.3};
    delta = 0.0001; tev = 0.001; tbeta = 1.0;
    output = start(p, delta, tev, tbeta);

    correct = {1.9662112, 2.06616, 0.8135182, 0.632185, 0.5964, 0.649624};
    correct_lambda_s = 41.517752; correct_t_eff = 1.0118507;

    THEN( "T1, debye waller, and effective temp are correctly computed" ){
      check( p, correct, output, correct_lambda_s, correct_t_eff );
    } // THEN

  } // GIVEN
  /*

  GIVEN( "more frequency distribution values" ){

    p = {.001, .003, .006, .008, .01, .03, .06, .08, 1.0, 3.0, 6.0, 8.0};
    tbeta = 2.0;

    WHEN( "spacing and temperature are very medium" ){
      delta = 0.015, tev = 1.723477E-2;
      output = start(p, delta, tev, tbeta);

      correct = {2.691905E-3, 4.031148E-3, 2.841214E-3, 2.247653E-3, 
                 2.014354E-3, 4.746877E-3, 7.851893E-3, 8.945384E-3, 
                 9.771149E-2, 0.2604205, 0.4686490, 0.5680046};
      correct_lambda_s = 0.238142827, correct_t_eff = 4.298571365;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        check( p, correct, output, correct_lambda_s, correct_t_eff );
      } // THEN
    } // WHEN

    WHEN( "temperature is changed" ){
      delta = 0.015; tev = 4.3086925E-2;
      output = start(p, delta, tev, tbeta);
      correct = {1.54443108E-2, 1.82883208E-2, 1.07199712E-2, 7.37431363E-3, 
                 5.96172958E-3, 1.30407541E-2, 2.04553735E-2, 2.24449206E-2, 
                 0.238766209, 0.62462916, 1.10947058, 1.332378};
      correct_lambda_s = 0.6485580; correct_t_eff = 1.836531025;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        check( p, correct, output, correct_lambda_s, correct_t_eff );
      } // THEN
    } // WHEN

    WHEN( "spacing is very large" ){
      tev = 4.3086925E-2; delta = 4.015;
      output = start(p, delta, tev, tbeta);

      correct = {2.36543334e-7, 2.20419880e-5, 2.20419880e-5, 1.95928782e-5, 
                 1.83683233e-5, 4.40839761e-5, 7.34732935e-5, 8.39694783e-5, 
                 9.18416169e-4, 2.44910978e-3, 4.40839761e-3, 5.34351226e-3};
      correct_lambda_s = 2.2081257E-3; correct_t_eff = 459.94244303;

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        check( p, correct, output, correct_lambda_s, correct_t_eff );
      } // THEN
    } // WHEN
  } // GIVEN
  */
} // TEST_CASE
