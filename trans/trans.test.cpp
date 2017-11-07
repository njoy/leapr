#define CATCH_CONFIG_MAIN
#include <iostream>
#include <vector>
#include "trans.h"
#include "../catch.hpp"

void equal( double a, double b ){
    if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
    if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

TEST_CASE( "trans" ){
    GIVEN( "inputs" ){
        std::vector<double> alpha {0.10,  0.20,  0.40, 0.50};
        std::vector<double> beta  {0.15,  0.18,  0.22};
        std::vector<double> temps {200.0}; 
        std::vector<std::vector<std::vector<double>>> sym_sab ( alpha.size(),
          std::vector<std::vector<double>> ( beta.size(),
            std::vector<double> ( temps.size(), 0.0 ) ) );
        sym_sab = { { { 0.1 }, { 0.2 }, { 0.3 } },
                    { { 0.4 }, { 0.5 }, { 0.6 } },
                    { { 0.7 }, { 0.8 }, { 0.9 } },
                    { { 1.0 }, { 1.1 }, { 1.2 } } };
        int lat = 3;
        double trans_weight = 0.03;
        double delta = 220.0;
        double diffusion_const = 1.5;
        double sc = 1.0;
        double arat = 1.0;
        double tev = 0.017;
        int itemp = 0;
        double lambda_s = 0.002;

        trans( alpha, beta, lat, trans_weight, delta, diffusion_const,
                sc, arat, tev, sym_sab, itemp, lambda_s );
        equal( sym_sab[0][0][0], 0.23049978 );
        equal( sym_sab[0][1][0], 0.25982880 );
        equal( sym_sab[1][0][0], 0.62197701 );
        equal( sym_sab[2][0][0], 1.08210491 );
        equal( sym_sab[3][0][0], 1.41011128 );


    } // GIVEN
} // TEST CASE

 //       sym_sab = { { { 0.1, 0.3 }, { 0.2, 0.4 }, { 0.3, 0.5 } },
 //                   { { 0.4, 0.6 }, { 0.5, 0.7 }, { 0.6, 0.8 } },
 //                   { { 0.7, 0.9 }, { 0.8, 1.0 }, { 0.9, 1.1 } },
 //                   { { 1.0, 1.2 }, { 1.1, 1.3 }, { 1.2, 1.4 } } };
 
