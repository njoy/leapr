#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "bfill.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}



TEST_CASE( "bfill" ){
  GIVEN( "inputs" ){
    WHEN( "the first beta value is tiny (<=1e-9)" ){
      int maxbb = 11;
      std::vector<double> bex   ( maxbb, 0.0 );
      std::vector<double> rdbex ( maxbb, 0.0 );
      std::vector<double> beta { 0.0, 0.15, 0.30, 0.60, 1.20 };
      std::vector<double> correct_bex {-1.2, -0.6, -0.3, -0.15, 0.0, 0.15, 0.3,
        0.6, 1.2, 0.0, 0.0};
      std::vector<double> correct_rdbex {1.6666667, 3.3333333, 6.6666667, 
        6.6666667, 6.6666667, 6.6666667, 3.3333333, 1.6666667, 0.0, 0.0, 0.0};
      THEN( "the vectors returned are correct" ){
        REQUIRE( bfill(bex, rdbex, beta) == 9 );
        equal_vec( bex, correct_bex );
        equal_vec( rdbex, correct_rdbex );
      } // THEN

      std::fill(bex.begin(), bex.end(), 0.0);
      std::fill(rdbex.begin(), rdbex.end(), 0.0);
      beta[0] = 1e-9;

      THEN( "the vectors returned are correct" ){
        REQUIRE( bfill(bex, rdbex, beta) == 9 );
        equal_vec( bex, correct_bex );
        equal_vec( rdbex, correct_rdbex );
      } // THEN
    } // WHEN
    WHEN( "the first beta value is of reasonable (>1e-9) size" ){
      int maxbb = 11;
      std::vector<double> bex   ( maxbb, 0.0 );
      std::vector<double> rdbex ( maxbb, 0.0 );
      std::vector<double> beta { 1.1e-9, 0.15, 0.30, 0.60, 1.20 };
      REQUIRE( bfill(bex, rdbex, beta) == 10 );
      std::vector<double> correct_bex {-1.2, -0.6, -0.3, -0.15, 0.0, 0.15, 
        0.3, 0.6, 1.2, 0.0, 0.0};
      std::vector<double> correct_rdbex {1.6666666, 3.3333333, 6.666666, 
        6.666666, 6.666666, 6.666666, 3.3333333, 1.6666666, 0.0, 0.0, 0.0};
      THEN( "the vectors returned are correct" ){
      //  equal_vec( bex, correct_bex );
      //  equal_vec( rdbex, correct_rdbex );
      } // THEN
      
    } // WHEN
  } // GIVEN
} // TEST CASE
