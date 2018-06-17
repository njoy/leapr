#include "catch.hpp"
/*
#include "iel/iel_util/terpa.h"


TEST_CASE( "terpa" ){
  double x, y, xnext;
  int idis, ip, ir;
  std::tuple<double,int> out;

  GIVEN( "small values in input vector" ){
    THEN( "110 multiple times then 150, defaulting to xbig" ){
      std::vector<double> a { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13, 21.21, 34.34 };
      x = 1; y = 2; ip = 2; ir = 1;
      out = terpa( y, x, a, ip, ir );
      REQUIRE( 1.0e12 == Approx( std::get<0>(out) ).epsilon(1e-6) );
    } // THEN

  } // GIVEN

  GIVEN( "x value is not equal to a value in the a vector" ){
    THEN( "110 then 120 and 130" ){
      std::vector<double> a (100);
      for ( int i = 0; i < 100; ++i ){ a[i] = 20 + i * 0.1; }
      x = 2; y = 1; ip = 2; ir = 1;
      out = terpa( y, x, a, ip, ir );
      REQUIRE( 24.6 == Approx( std::get<0>(out) ).epsilon(1e-6) );


    } // THEN

  } // GIVEN
  GIVEN( "x is equal to a value in the a vector " ){
    THEN( "110 then 120 and 140" ){
      std::vector<double> a (100);
      for ( int i = 0; i < 100; ++i ){ a[i] = 20 + i * 0.1; }
      x = 25; y = 1; ip = 2; ir = 1;
      out = terpa( y, x, a, ip, ir );
      REQUIRE( 25.2 == Approx( std::get<0>(out) ).epsilon(1e-6) );

    } // THEN

  } // GIVEN
  GIVEN( " " ){
    THEN( "110 immediately to 150" ){
      std::vector<double> a { 2.01, 2.02, 2.03, 2.05, 2.08, 2.13, 2.21, 2.34, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      x = 12.5; y = 1; ip = 2; ir = 1;
      out = terpa( y, x, a, ip, ir );
      REQUIRE( 1.0e12 == Approx( std::get<0>(out) ).epsilon(1e-6) );
      std::cout << std::get<0>(out) << std::endl;

    } // THEN

  } // GIVEN

} // TEST CASE
*/
