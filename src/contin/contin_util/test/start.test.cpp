#include "catch.hpp"
#include "contin/contin_util/start.h"
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>

void equal2( double a, double b ){
    REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal2_vec( std::vector<double> a, std::vector<double> b ){
    REQUIRE( a.size() == b.size() );
    for ( int i = 0; i < a.size(); ++i ){
        equal2( a[i], b[i] );
    }
}

TEST_CASE( "start function" ){
  GIVEN( "phonon distribution, spacing, and temp in ev" ){
    std::vector<double> p {0.1, 0.2, 0.3, 0.5, 0.8, 1.3};
    double delta = 0.0001, tev = 0.001, tbeta = 1.0;
    std::vector<double> correct {1.96621124173, 2.06616003278, 
      0.813518278800, 0.632185383948, 0.596399990773, 0.649624806321};
    double correct_lambda_s = 41.517752;
    double correct_t_eff    = 1.0118507;
    THEN( "all entries of output have at most 1e-6 percentage error" ){
      std::tuple<double,double> output = start(p, delta, tev, tbeta);
      equal2( std::get<0>(output), correct_lambda_s );
      equal2( std::get<1>(output), correct_t_eff );
      equal2_vec( p, correct );
    } // THEN
    p = {.001, .003, .006, .008, .01, .03, .06, .08, 1.0, 3.0, 6.0, 8.0};
    delta = 0.015, tev = 1.723477E-2, tbeta = 2.0;

    correct = {2.69190582E-3, 4.03114887E-3, 2.84121466E-3, 
      2.24765366E-3, 2.01435442E-3, 4.74687728E-3, 7.85189344E-3, 
      8.94538406E-3, 9.77114977E-2, 0.260420592, 0.468649045, 0.568004629};
    correct_lambda_s = 0.238142827;
    correct_t_eff = 4.298571365;

    THEN( "all entries of output have at most 1e-6 percentage error" ){
      std::tuple<double,double> output = start(p, delta, tev, tbeta);
      equal2( std::get<0>(output), correct_lambda_s );
      equal2( std::get<1>(output), correct_t_eff );
      equal2_vec( p, correct );
    } // THEN

    tev = 4.3086925E-2;

    correct = {1.54443108E-2, 1.82883208E-2, 1.07199712E-2, 7.37431363E-3, 
               5.96172958E-3, 1.30407541E-2, 2.04553735E-2, 2.24449206E-2, 
               0.238766209, 0.62462916, 1.10947058, 1.332378};
    correct_lambda_s = 0.6485580;
    correct_t_eff = 1.836531025;

    THEN( "all entries of output have at most 1e-6 percentage error" ){
      std::tuple<double,double> output = start(p, delta, tev, tbeta);
      equal2( std::get<0>(output), correct_lambda_s );
      equal2( std::get<1>(output), correct_t_eff );
      equal2_vec( p, correct );
    } // THEN

    delta = 4.015;

    correct = {2.36543334e-7, 2.20419880e-5, 2.20419880e-5, 1.95928782e-5, 
               1.83683233e-5, 4.40839761e-5, 7.34732935e-5, 8.39694783e-5, 
               9.18416169e-4, 2.44910978e-3, 4.40839761e-3, 5.34351226e-3};
    correct_lambda_s = 2.2081257E-3;
    correct_t_eff = 459.94244303;

    THEN( "all entries of output have at most 1e-6 percentage error" ){
      std::tuple<double,double> output = start(p, delta, tev, tbeta);
      equal2( std::get<0>(output), correct_lambda_s );
      equal2( std::get<1>(output), correct_t_eff );
      equal2_vec( p, correct );
    } // THEN
  } // GIVEN
} // TEST_CASE
