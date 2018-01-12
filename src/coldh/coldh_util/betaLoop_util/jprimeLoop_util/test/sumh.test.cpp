#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../sumh.h"

void equal( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}



TEST_CASE( "sumh" ){
  GIVEN( "inputs" ){
    WHEN( "j == 0" ){
      equal( sumh(0,1,5e-3), 2.77776376E-6 );
      equal( sumh(0,3,5e-3), 1.41722943E-18 );
      equal( sumh(0,8,1e-5), 8.42138983E-96 );
    }
    WHEN( "neither j or jp equal 0" ){
      equal( sumh(1,3,1.5), 6.96384680E-3 );
      equal( sumh(1,1,5e-3), 0.33333055 );
      equal( sumh(3,8,1e-5), 2.3450205E-59 );
    }
    WHEN( "j doesnt equal 0 and jp == 0" ){
      equal( sumh(1,0,1.5),  0.156953021 );
      equal( sumh(1,0,5e-3), 2.77776376E-6 );
      equal( sumh(1,0,0.35), 1.32811188E-2 );
    }



  } // THEN
} // TEST CASE
