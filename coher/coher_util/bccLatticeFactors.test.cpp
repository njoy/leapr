#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "bccLatticeFactors.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( abs(b-a) < 1e-6 );
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


TEST_CASE( "Function to Compute BCC Lattice Factors" ){
  GIVEN( "inputs" ){
    THEN( "outputs" ){
      REQUIRE( true );

    } // THEN
  } // GIVEN
} // TEST CASE

