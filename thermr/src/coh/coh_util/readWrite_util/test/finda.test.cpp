#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coh/coh_util/readWrite_util/finda.h"


TEST_CASE( "finda" ){
  std::fstream file("temp_float");
  int i, na, nbuf = 10;
  std::vector<double> buf(nbuf,0.0);

  GIVEN( "inow value of 0" ){
    i = 1;
    WHEN( "a vector is small" ){
      na = 2;
      std::vector<double> a(na,0.0);

      finda( i, na, file, a, buf, nbuf );

      THEN( "values of a are read in from file correctly" ){
        REQUIRE( 1.0 == Approx(a[0]).epsilon(1e-6) );
        REQUIRE( 2.0 == Approx(a[1]).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "a vector is medium" ){
      na = 5;
      std::vector<double> a(na,0.0);

      finda( i, na, file, a, buf, nbuf );

      THEN( "values of a are read in from file correctly" ){

        REQUIRE( 1.0 == Approx(a[0]).epsilon(1e-6) );
        REQUIRE( 2.0 == Approx(a[1]).epsilon(1e-6) );
        REQUIRE( 3.0 == Approx(a[2]).epsilon(1e-6) );
        REQUIRE( 4.0 == Approx(a[3]).epsilon(1e-6) );
        REQUIRE( 5.0 == Approx(a[4]).epsilon(1e-6) );

      } // THEN

    } // WHEN
  } // GIVEN

  GIVEN( "inow doesn't equal 1" ){
    i = 2;

    WHEN( "buf is initialized to all zeros originally" ){
      na = 5;
      std::vector<double> a(na,0.0);

      finda( i, na, file, a, buf, nbuf );

      THEN( "buf is not read in" ){
        for ( auto& entry : a ){
          REQUIRE( 0.0 == Approx(entry).epsilon(1e-6) );
        }
      } // THEN

    } // WHEN

    WHEN( "buf is populated with other values" ){
      buf = { 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.0, 10.1 };
      THEN( "buf is not read in, a is filled with existent buf values" ){
        { 
        na = 2;
        std::vector<double> a(na,0.0);

        finda( i, na, file, a, buf, nbuf );

        REQUIRE( 3.4 == Approx(a[0]).epsilon(1e-6) );
        REQUIRE( 4.5 == Approx(a[1]).epsilon(1e-6) );
        }
        { 
        na = 5;
        std::vector<double> a(na,0.0);

        finda( i, na, file, a, buf, nbuf );

        REQUIRE( 6.7  == Approx(a[0]).epsilon(1e-6) );
        REQUIRE( 7.8  == Approx(a[1]).epsilon(1e-6) );
        REQUIRE( 8.9  == Approx(a[2]).epsilon(1e-6) );
        REQUIRE( 9.0  == Approx(a[3]).epsilon(1e-6) );
        REQUIRE( 10.1 == Approx(a[4]).epsilon(1e-6) );
        }

      } // THEN

    } // WHEN

  } // GIVEN

} // TEST CASE
