#include "catch.hpp"
#include "discre/discre_util/bfill.h"
#include "generalTools/testing.h"


TEST_CASE( "bfill" ){
  GIVEN( "inputs" ){
    WHEN( "the first beta value is tiny (<=1e-9)" ){
      int maxbb = 11;
      std::vector<double> 
        rdbex (maxbb, 0.0),
        beta {0.0,0.15,0.3,0.6,1.2},
        correct_bex {-1.2,-0.6,-0.3,-0.15,0.0,0.15,0.3,0.6,1.2,0.0,0.0},
        correct_rdbex {1.6666667,3.3333333,6.6666667,6.6666667,6.6666667,
                       6.6666667,3.3333333,1.6666667,0.0,0.0,0.0};
      THEN( "the vectors returned are correct" ){
        auto output = bfill(maxbb, rdbex, beta);
        REQUIRE( std::get<0>(output) == 9 );
        REQUIRE( ranges::equal(std::get<1>(output),correct_bex,equal) );
        REQUIRE( ranges::equal(rdbex,correct_rdbex,equal) );
      } // THEN

      std::fill(rdbex.begin(), rdbex.end(), 0.0);
      beta[0] = 1e-9;

      THEN( "the vectors returned are correct" ){
        auto output = bfill(maxbb, rdbex, beta);
        REQUIRE( std::get<0>(output) == 9 );
        REQUIRE( ranges::equal(std::get<1>(output),correct_bex,equal) );
        REQUIRE( ranges::equal(rdbex,correct_rdbex,equal) );
      } // THEN
    } // WHEN
    WHEN( "the first beta value is of reasonable (>1e-9) size" ){
      int maxbb = 11;
      std::vector<double> 
        rdbex (maxbb, 0.0),
        beta {1.1e-9,0.15,0.30,0.60,1.20},
        correct_bex {-1.2,-0.6,-0.3,-0.15,-1.1e-9,1.1e-9,0.15,0.3,0.6,1.2,0.0},
        correct_rdbex {1.6666666,3.3333333,6.666666,6.6666666,4.545455e8, 
                       6.6666666,6.66666666,3.3333333,1.6666666,0.0,0.0 };
      auto output = bfill(maxbb, rdbex, beta);
      REQUIRE( std::get<0>(output) == 10 );
      THEN( "the vectors returned are correct" ){
        REQUIRE( ranges::equal(std::get<1>(output),correct_bex,equal) );
        REQUIRE( ranges::equal(rdbex,correct_rdbex,equal) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE


