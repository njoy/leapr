#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../catch.hpp"
#include "coher.h"

void equal( double a, double b ){
  std::cout << a << "    " << b << std::endl;

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



TEST_CASE( "coher" ){
  REQUIRE( true );
} // TEST CASE
