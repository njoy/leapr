#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "hexLatticeFactorsMid.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( std::abs(b-a) < 1e-6 );
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


TEST_CASE( "Function to Compute Hex Lattice Factors" ){
  double a = 1e-9, c1 = 1.5e15, c2 = 2.5e15, tsqx = 9.6e17,
    t2 = 3.5e-5, ulim = 9.6e19, c = 3.58e-8, tsq = 0, wint = 0;
  int i = 0, ifl = 1, lat = 3, nw = 60000, imax = 5;
  std::vector<double> b (60000, 0.0);

  int i1 = 1, l1 = 0, i2m = 2; 
  GIVEN( "few iterations" ){
    THEN( "outputs" ){
      hexLatticeFactorsMid( a, tsq, c1, c2, lat, nw, tsqx, b, ifl, 
          i, wint, t2, ulim, imax, c, i1, i2m, l1 );
      for ( auto i = 0; i < 30; ++i ){ std::cout << b[i] << std::endl; }

    } // THEN
  } // GIVEN
} // TEST CASE

