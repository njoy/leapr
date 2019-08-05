#include "catch.hpp"
#include "coldh/coldh_util/terpk.h"
#include "generalTools/interpolate.h"
#include <range/v3/all.hpp>


TEST_CASE( "terpk" ){
  std::vector<double> ska(6);
  double delta, be1, be2, be3;

  GIVEN( "desired value within the range of the ska vector" ){

    WHEN( "desired value is not a point in ska vector" ){
      ska = {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
      delta = 0.2; be1 = 0.6; be2 = 0.01;
      THEN( "interpolation is correctly performed" ){
        REQUIRE( terpk( ska, delta, be1 ) == Approx(5.5).epsilon(1e-6) );
        REQUIRE( terpk( ska, delta, be2 ) == Approx(1.155).epsilon(1e-6) );
        //auto xVector = ranges::view::iota(0,int(ska.size())) | ranges::view::transform([delta](auto x){return x*delta;});;
        //std::cout << xVector << std::endl;
        //std::cout << interpolate( ranges::view::zip(xVector,ska), be1 ) << std::endl;
      } // THEN
    } // WHEN

    WHEN( "desired value is a point in ska vector" ){
      ska = {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
      delta = 0.2; be1 = 0.0; be2 = 0.2; be3 = 0.8;
      THEN( "corresponding point is returned" ){
        REQUIRE( terpk( ska, delta, be1 ) == Approx(1.1).epsilon(1e-6) );
        REQUIRE( terpk( ska, delta, be2 ) == Approx(2.2).epsilon(1e-6) );
        REQUIRE( terpk( ska, delta, be3 ) == Approx(8.8).epsilon(1e-6) );
      } // THEN
    } // WHEN

  } // GIVEN

  GIVEN( "desired value is greater than the range of the ska vector" ){
    ska = {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
    delta = 0.2; be1 = 1.0; be2 = 5.01;
    THEN( "a value of 1.0 is returned" ){
      REQUIRE( terpk( ska, delta, be1 ) == Approx(1.0).epsilon(1e-6) );
      REQUIRE( terpk( ska, delta, be2 ) == Approx(1.0).epsilon(1e-6) );
    } // THEN
  } // GIVEN
} // TEST CASE
