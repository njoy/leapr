#define CATCH_CONFIG_MAIN
#include "catch.hpp"
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

                std::vector<double> output = convol( t1, t2, delta, t2.size(), t1.size() );
                std::vector<double> correct = {5.499735535, 3.341519474, 
                  1.243674184, 0.6360538183, 0.4786883192, 0.4560000019, 
                  0.6960000001, 1.103999981, 2.039999954, 2.879999935, 
                  1.919999957, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                equal_vec(output, correct);
           
            }
            {
                std::vector<double> t1 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21};
                std::vector<double> t2 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21,
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                                          0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
                double delta = 0.5;
                std::vector<double> output = convol( t1, t2, delta, t2.size(), t1.size() );
                std::vector<double> correct {1.383467893E-2, 1.489987296E-2, 
                  1.52705054E-2, 1.60847176E-2, 1.59448573E-2, 1.94499999E-2, 
                  2.88499997E-2, 3.64999997E-2, 3.58999986E-2, 3.35999982E-2, 
                  2.204999862E-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                equal_vec(output, correct);
            }
            {
                std::vector<double> t1 = {0.001,0.003,0.006,0.01,0.04,0.06,0.08,
                 0.096,0.11,0.15};
                std::vector<double> t2 = {0.001,0.003,0.006,0.01,0.04,0.06,0.08,
                 0.096,0.11,0.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                double delta = 0.2;
                std::vector<double> output = convol( t1, t2, delta, t2.size(), t2.size() );
                std::vector<double> correct = {4.99725082E-3, 4.54868791E-3, 
                  4.11462771E-3, 3.43174809E-3, 2.55677869E-3, 1.66090654E-3, 
                  8.36569270E-4, 8.06987691E-4, 1.08057153E-3, 1.73239998E-3, 
                  2.82799995E-3, 4.25599996E-3, 5.94399996E-3, 8.11200004E-3, 
                  8.96320013E-3, 9.02400023E-3, 8.18000037E-3, 6.60000032E-3, 
                  4.50000042E-3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                  0.0, 0.0};
                equal_vec(output, correct);

            }
        } // THEN
    } // GIVEN
} // TEST CASE

