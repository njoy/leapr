#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "hexLatticeFactorsHelper.h"

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


TEST_CASE( "Function to Compute Hexagonal Lattice Factors" ){
  int lat = 2, k = 0, ifl = 1, i = 0, nw = 6;
  double tsq = 0.1, tsqx  = 9.6, wint = 0, eps = 5e-5, f = 3;
  std::vector<double> b ( 60000, 0.0 );

  GIVEN( "k is not positive or tsq <= tsqx" ){
    THEN( "b's first two entries are populted with tsq and f, respectively" ){
      hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, i, wint, nw, eps, f );

      equal( b[0], 0.1 ); equal( b[1], 3.0 ); equal( k, 1 ); equal( i, 0 );
      for ( auto i = 2; i < b.size(); ++i ){ equal( b[i], 0 ); }

      f = 100; k = 0;

      hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, i, wint, nw, eps, f );

      equal( b[0], 0.1 ); equal( b[1], 100 ); equal( k, 1 ); equal( i, 0 );
      for ( auto i = 2; i < b.size(); ++i ){ equal( b[i], 0 ); }


      tsq = 50; k = 0;

      hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, i, wint, nw, eps, f );
      equal( b[0], 50 ); equal( b[1], 100 ); equal( k, 1 ); equal( i, 0 );
      for ( auto i = 2; i < b.size(); ++i ){ equal( b[i], 0 ); }


    } // THEN
  } // GIVEN

  GIVEN( "k is positive and tsq > tsqx" ){
    WHEN( "loop is able to finish" ){
      THEN( "we're in the first situation" ){
        tsq = 30.1, tsqx = 25;
        k = 1; f = 4.5e-2;
        hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, i, wint, nw, eps, f );

//        equal( b[0], 0 ); equal( b[1], 0 ); equal( b[2], 700 ); 
//        equal( b[3], 4.5e-2); equal( k, 2 ); equal( i, 1 );
//        for ( auto i = 4; i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
    
    WHEN( "loop isn't able to finish" ){
      THEN( "we're in the second situation" ){
        tsq = 700, tsqx = 500;
        k = 1; f = 4.5e-2;
        hexLatticeFactorsHelper( k, tsq, tsqx, b, ifl, i, wint, nw, eps, f );

        equal( b[0], 0 ); equal( b[1], 0 ); equal( b[2], 700 ); 
        equal( b[3], 4.5e-2); equal( k, 2 ); equal( i, 1 );
        for ( auto i = 4; i < b.size(); ++i ){ equal( b[i], 0 ); }


      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE

