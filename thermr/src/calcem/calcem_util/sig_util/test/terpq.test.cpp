#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calem_util/sig_util/terpq.h"


TEST_CASE( "terpq" ){
  double out;
  GIVEN( "desired x value is below provided range (x < x1)" ){
    WHEN( "x < x1 and y1 <= y2" ){
      THEN( "terp1 is called, interpolation code 3 is used" ){
        //                x1    x2    x3         x       y1     y2     y3   
        double in[8] = { 2.25, 2.30, 2.35, 5.0428643E-2, 0.10, 1.40, 2.10 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( -224.55032636 == Approx( out ).epsilon(1e-6) );

        in[3] = 2.24;
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( -0.1634642073 == Approx( out ).epsilon(1e-6) );

        in[4] = 1.10; in[5] = 2.40; in[6] = 9.80;
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 0.83653579541 == Approx( out ).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "x < x1 and y1 > y2" ){
      THEN( "y = y1 is returned" ){
        //                x1    x2    x3         x       y1     y2     y3    y
        double in[8] = { 2.25, 2.30, 2.35, 5.0428643E-2, 4.10, 3.40, 2.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 4.1 == Approx( out ).epsilon(1e-6) );

        in[3] = 2.24; in[4] = 9.10; in[5] = 8.40; in[6] = 2.80;
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 9.1 == Approx( out ).epsilon(1e-6) );

      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( "desired x value is above provided range (x > x3)" ){

    WHEN( "x > x3 and y3 <= y2" ){
      THEN( "terp1 is called, interpolation code 2 is used" ){
        //                x1    x2    x3    x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 15.0, 4.10, 3.40, 2.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 1.155395387 == Approx( out ).epsilon(1e-6) );

        in[3] = 9.96;
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 2.098129357 == Approx( out ).epsilon(1e-6) );

        in[3] = 9.99; in[4] = 1.30; in[5] = 1.20; in[6] = 1.00;
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 0.998848921 == Approx( out ).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "x > x3 and y3 > y2" ){
      THEN( "y is set to y3 and returned" ){
        //                x1    x2    x3    x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 15.0, 1.10, 3.40, 8.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 8.10 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( "desired x value is within range (x1 < x < x3)" ){
    WHEN( "big steps are used" ){
      AND_WHEN( "x >= x2" ){
        THEN( "terp1 is used, interpolating between x2 and x3" ){
          //                x1    x2    x3   x    y1    y2     y3    y
          double in[8] = { 1.05, 3.00, 9.95, 5.0, 1.10, 3.40, 8.10, 0.0 };
          out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
          REQUIRE( 4.7525182 == Approx( out ).epsilon(1e-6) );
        } // THEN
      } // WHEN
      AND_WHEN( "big steps and x < x2" ){
        THEN( "terp1 is used, interpolating between x1 and x2" ){
          //                x1    x2    x3   x    y1    y2     y3    y
          double in[8] = { 1.05, 3.00, 9.95, 2.0, 1.10, 3.40, 8.10, 0.0 };
          out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
          REQUIRE( 2.2205129 == Approx( out ).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "small steps" ){
      THEN( "terp1 is used, interpolating between x1 and x2" ){
        //                x1    x2    x3   x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 2.0, 1.10, 1.42, 1.73, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6] );
        REQUIRE( 1.268652857 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
 
