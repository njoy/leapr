#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "generalTools/trapezoidIntegral.h"
#include "range/v3/all.hpp"


TEST_CASE( "trapezoid integral" ){
  GIVEN( "Zipped XY input range and a func to act on X and Y" ){ 
    using std::pow;
    auto func1 = [](auto xy ){ 
      return pow(std::get<0>(xy),2.0) + std::get<1>(xy)*2.5 - 1; };
    auto func2 = [](auto xy){ 
      return pow(std::get<0>(xy),3.5) + pow(std::get<1>(xy),1.2); };
    auto func3 = [](auto xy){ 
      return pow(std::get<0>(xy),3.5) * pow(std::get<1>(xy),0.02); };

    WHEN( "Input x grid is uniform" ){
      std::vector<double> x {0, 1, 2, 3, 4},
                          y {0, 4, 3, 1, 2};
      auto XY = ranges::view::zip(x,y);
      THEN( "func is applied and trapezoid integral is taken" ){
        REQUIRE(trapezoidIntegral(XY,func1) == Approx(40.5).epsilon(1e-6));
        REQUIRE(trapezoidIntegral(XY,func2) == Approx(134.243).epsilon(1e-6));
        REQUIRE(trapezoidIntegral(XY,func3) == Approx(124.25194).epsilon(1e-6));
      } // THEN
    } // WHEN

    WHEN( "Input x grid is non uniform" ){
      std::vector<double> x {0, 1, 2, 4, 8, 16, 32},
                          y {1, 3, 7, 3, 2,  4, 5 };
      auto XY = ranges::view::zip(x,y);
      THEN( "func is applied and trapezoid integral is taken" ){
        REQUIRE(trapezoidIntegral(XY,func1) == Approx(11978.5).epsilon(1e-6));
        REQUIRE(trapezoidIntegral(XY,func2) == Approx(1688772.5).epsilon(1e-6));
        REQUIRE(trapezoidIntegral(XY,func3) == Approx(1742776.8).epsilon(1e-6));
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
