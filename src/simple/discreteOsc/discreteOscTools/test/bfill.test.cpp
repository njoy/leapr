#include "catch.hpp"
#include "simple/discreteOsc/discreteOscTools/bfill.h"
#include "simple/generalTools/testing.h"

TEST_CASE( "bfill" ){
  GIVEN( "inputs" ){
    WHEN( "the first beta value is tiny (<=1e-9)" ){
      int maxbb = 11;
      std::vector<double> bex (maxbb, 0.0), rdbex (maxbb, 0.0),
        beta { 0.0, 0.15, 0.30, 0.60, 1.20 },
        correct_bex {-1.2,-0.6,-0.3,-0.15,0.0,0.15,0.3,0.6,1.2,0.0,0.0},
      correct_rdbex {1.6666667, 3.3333333, 6.6666667, 
        6.6666667, 6.6666667, 6.6666667, 3.3333333, 1.6666667, 0.0, 0.0, 0.0};

      THEN( "the vectors returned are correct" ){
        REQUIRE( bfill(bex, rdbex, beta) == 9 );
        checkVec( bex, correct_bex );
        checkVec( rdbex, correct_rdbex );
      } // THEN

      std::fill(bex.begin(), bex.end(), 0.0);
      std::fill(rdbex.begin(), rdbex.end(), 0.0);
      beta[0] = 1e-9;

      THEN( "the vectors returned are correct" ){
        REQUIRE( bfill(bex, rdbex, beta) == 9 );
        checkVec( bex, correct_bex );
        checkVec( rdbex, correct_rdbex );
      } // THEN
    } // WHEN
    WHEN( "the first beta value is of reasonable (>1e-9) size" ){
      int maxbb = 11;
      std::vector<double> bex (maxbb, 0.0), rdbex (maxbb, 0.0), 
        beta { 1.1e-9, 0.15, 0.30, 0.60, 1.20 },
        correct_bex {-1.2,-0.6,-0.3,-0.15,-1.1e-9,1.1e-9,0.15,0.3,0.6,1.2,0.0},
        correct_rdbex {1.6666666, 3.3333333, 6.666666, 6.6666666, 4.545455e8, 
          6.6666666, 6.66666666, 3.3333333, 1.6666666, 0.0, 0.0 };

      REQUIRE( bfill(bex, rdbex, beta) == 10 );
      THEN( "the vectors returned are correct" ){
        checkVec( bex, correct_bex );
        checkVec( rdbex, correct_rdbex );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
