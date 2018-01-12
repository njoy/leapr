#include "catch.hpp"
#include "coldh/coldh_util/terpk.h"

void equal1( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "terpk" ){
  GIVEN( "desired value within the range of the ska vector" ){
    WHEN( "desired value is not a point in ska vector" ){
      std::vector<double> ska {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
      double delta = 0.2;
      double be1 = 0.6;
      double be2 = 0.01;
      THEN( "interpolation is correctly performed" ){
        equal1( terpk( ska, delta, be1 ), 5.5 );
        equal1( terpk( ska, delta, be2 ), 1.155 );
      } // THEN
    } // WHEN
    WHEN( "desired value is a point in ska vector" ){
      std::vector<double> ska {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
      double delta = 0.2;
      double be1 = 0.0;
      double be2 = 0.2;
      double be3 = 0.8;
      THEN( "corresponding point is returned" ){
        equal1( terpk( ska, delta, be1 ), 1.1 );
        equal1( terpk( ska, delta, be2 ), 2.2 );
        equal1( terpk( ska, delta, be3 ), 8.8 );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "desired value is greater than the range of the ska vector" ){
    std::vector<double> ska {1.1, 2.2, 3.3, 5.5, 8.8, 13.13}; 
    double delta = 0.2;
    double be1 = 1.0;
    double be2 = 5.01;
    THEN( "a value of 1.0 is returned" ){
      equal1( terpk( ska, delta, be1 ), 1.0 );
      equal1( terpk( ska, delta, be2 ), 1.0 );
    } // THEN
  } // GIVEN
} // TEST CASE
