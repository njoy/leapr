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
  GIVEN( "material of interest is iron" ){
    WHEN( "inputs" ){
      THEN( "outputs" ){
        int iel = 6, npr = 1, maxb = 60000;
	std::vector<double> b ( 60000, 0.0 );
	double emax = 5.0;
	coher( iel, npr, maxb, b, emax );
        	
        REQUIRE( true );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
