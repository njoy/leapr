#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "hexagonalFactorsHelper.h"

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
  int lat = 2, l1 = 0, l2 = 0, l3 = 0, k = 0, ifl = 1, i = 0, nw = 6;
  double w = 0, tsq = 0.1, tsqx  = 9.6, ulim = 9.6e19, t2 = 3.5e-5, w1 = 1,
         w2 = 0.5, w3 = 1, wint = 0, eps = 5e-5;
  std::vector<double> b ( 60000, 0.0 );

  hexagonalLatticeFactorsHelper( lat, l1, l2, l3, w, k, tsq, tsqx, b, ifl, i,
      ulim, t2, w1, w2, w3, wint, nw, eps );

  equal( b[0], 0.1 );
  equal( b[1], 3.1622776 );
  equal( w,    1.5811388 );
  equal( k, 1 );
  equal( i, 0 );
  for ( auto i = 2; i < b.size(); ++i ){ equal( b[i], 0 ); }

} // TEST CASE




TEST_CASE( "coher" ){
  REQUIRE( true );
} // TEST CASE
