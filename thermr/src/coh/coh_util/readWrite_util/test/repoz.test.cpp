#define CATCH_CONFIG_MAIN
#include "coh/coh_util/readWrite_util/repoz.h"


TEST_CASE( "repoz" ){
  GIVEN( "" ){
    std::fstream infile;
    infile.open( "temp_int" );
    std::vector<int> a(9);
    WHEN( "we read in the first part of the file" ){ 
      infile >> a[0] >> a[1] >> a[2];
      infile >> a[3] >> a[4] >> a[5];

      AND_WHEN( "we don't call repoz" ){
        THEN( "when we read, we pick up where we left off" ){
          infile >> a[6] >> a[7] >> a[8];
          REQUIRE( a[6] != a[0] );
          REQUIRE( a[7] != a[1] );
          REQUIRE( a[8] != a[2] );
        } // THEN
      } // AND WHEN

      AND_WHEN( "we call repoz" ){
        repoz( infile );
        THEN( "when we start reading again, we start at beginning" ){
          infile >> a[6] >> a[7] >> a[8];
          REQUIRE( a[6] == a[0] );
          REQUIRE( a[7] == a[1] );
          REQUIRE( a[8] == a[2] );
        } // THEN
      } // AND WHEN
    } // WHEN
    infile.close();
  } // GIVEN
} // TEST CASE
