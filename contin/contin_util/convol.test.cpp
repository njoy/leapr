#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "convol.h"
#include <iostream>

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

TEST_CASE( "convol" ){
    GIVEN( "two vectors" ){
        THEN( "correctly convolves the two" ){
            {
                std::vector<double> t1 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0};

                std::vector<double> t2 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                double delta = 0.03;

                std::vector<double> output = convol( t1, t2, delta );
                std::vector<double> correct = {3.8459762, 2.6993367, 1.0195307, 0.53364442, 0.37281623, 0.38400000, 0.62399999, 1.0079999, 1.7999999, 2.1599999, 0.95999997854232788, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                equal_vec(output, correct);
           
            }
            {
                std::vector<double> t1 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21};
                std::vector<double> t2 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21,
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
                double delta = 0.5;
                std::vector<double> output = convol( t1, t2, delta );
                std::vector<double> correct {1.1974704E-2, 1.3563056E-2, 1.3531928E-2, 1.3796487E-2, 1.3871143E-2, 1.7874999E-2, 2.6749999E-2, 3.1774999E-2, 3.0124998E-2, 2.5199998E-2, 1.1024999E-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                equal_vec(output, correct);
            }
        } // THEN
    } // GIVEN
} // TEST CASE

