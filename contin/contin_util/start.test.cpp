#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "start.h"
#include <vector>
#include <memory>
#include <iostream>

void equal( double a, double b ){
    REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
    REQUIRE( a.size() == b.size() );
    for ( int i = 0; i < a.size(); ++i ){
        equal( a[i], b[i] );
    }
}

TEST_CASE( "start function" ){
    GIVEN( "phonon distribution, spacing, and temp in ev" ){
        std::vector<double> p {0.1, 0.2, 0.3, 0.5, 0.8, 1.3 };
        std::vector<double> correct {1.96621124173, 2.06616003278, 0.813518278800, 0.632185383948, 0.596399990773, 0.649624806321};
        double correct_lambda_s = 41.51775265;
        double delta = 0.0001, tev = 0.001, tbeta = 1.0;
        THEN( "all entries of output have at most 1e-6 percentage error" ){
            equal( start(p, delta, tev, tbeta ), correct_lambda_s );
            equal_vec( p, correct );
        } // THEN
    } // GIVEN
} // TEST_CASE
