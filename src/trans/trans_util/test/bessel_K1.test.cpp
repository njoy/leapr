#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "trans/trans_util/bessel_K1.h"
#include <range/v3/all.hpp>


TEST_CASE( "bessel_K1 generating function" ){
  GIVEN( "values of x <= 1 at which to be evaluated" ){
        
      auto xRange = ranges::view::iota(0,20)
                  | ranges::view::transform([](int i){return 0.05*i + 5e-3;});
      std::vector<double> k1Vals = {199.985214, 18.085074, 9.3728868, 
        6.2586399, 4.6508874, 3.6659691, 2.9991255, 2.5169911, 2.1518758, 
        1.8657226, 1.6354620, 1.4462794, 1.2882126, 1.1543071, 1.0395568, 
        0.9402608, 0.8536217, 0.7774836, 0.7101568, 0.6502969 };

      THEN( "the modified bessel function value is returned" ){
        RANGES_FOR( auto t, ranges::view::zip(xRange,k1Vals) ){
      //    REQUIRE( std::get<1>(t) == 
      //             Approx(bessel_K1_gen(std::get<0>(t))).epsilon(1e-6) );
        }
      } // THEN
  } // GIVEN
  GIVEN( "values of x > 1 at which to be evaluated" ){

   /*
      auto xRange = ranges::view::iota(0,6)
                  | ranges::view::transform([](int i){return 5.0*i + 1.005;});
          std::vector<double> k1Vals = { 0.5968202, 0.0013366, 6.4867716e-6, 
          3.5885809e-8, 2.099515609e-10, 1.267212993e-12, 7.802144033e-15 };
	  RANGES_FOR( auto t, ranges::view::zip(xRange,k1Vals) ){
           std::cout << std::get<0>(t) << "    " << std::get<1>(t) << std::endl;
          REQUIRE( std::get<1>(t) == 
                   Approx(bessel_K1_gen(std::get<0>(t))).epsilon(1e-6) );
        } */


      std::vector<double> x_vec { 1.005, 1.5, 5.0, 10.0, 100.0 };
      std::vector<double> k1Vec{1.6304576457, 1.2431658735, 0.6002738602, 
                                     0.4107665699, 0.1257999504}; 
      THEN( "the modified bessel function value is returned, except the "
            "exponential piece is omitted" ){
      for( size_t i = 0; i < x_vec.size(); ++i ){
      //  REQUIRE(k1Vec[i] == Approx(bessel_K1_gen(x_vec[i]).epsilon(1e-6)));
      } // for
    } // THEN
  } // GIVEN
} // TEST CASE
