#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "trans/trans_util/bessel_K1.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

TEST_CASE( "bessel_K1 generating functin" ){
  GIVEN( "values of x <= 1 at which to be evaluated" ){
    std::vector<double> x_vec   { 0.005, 0.01, 0.1, 0.7, 0.999, 1.0 },
                        correct { 199.98521432, 99.9738941, 9.85384478, 
                                  1.05028353, 0.60293127, 0.6019072 }; 
    THEN( "the modified bessel function value is returned" ){
      auto outputs = ranges::view::iota(0,int(x_vec.size())) 
                   | ranges::view::transform( [&x_vec](int i){
                       return bessel_K1_gen(x_vec[i]); } );
      REQUIRE( ranges::equal(outputs,correct,equal) );
    } // THEN
  } // GIVEN
  GIVEN( "values of x > 1 at which to be evaluated" ){
    std::vector<double> x_vec { 1.005, 1.5, 5.0, 10.0, 100.0 },
                      correct {1.6304576457, 1.2431658735, 0.6002738602, 
                               0.4107665699, 0.1257999504}; 
    THEN( "the modified bessel function value is returned" ){
      auto outputs = ranges::view::iota(0,int(x_vec.size())) 
                   | ranges::view::transform( [&x_vec](int i){
                       return bessel_K1_gen(x_vec[i]); } );
      REQUIRE( ranges::equal(outputs,correct,equal) );
    } // THEN
  } // GIVEN
} // TEST CASE
