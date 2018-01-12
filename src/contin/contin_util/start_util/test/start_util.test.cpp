#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/start_util.h"
#include <vector>
#include <cmath>

void equal( double a, double b ){
    REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
    REQUIRE( a.size() == b.size() );
    for ( int i = 0; i < a.size(); ++i ){
        equal( a[i], b[i] );
    }
}

TEST_CASE( "fsum" ){
    WHEN( "n = 1 (used for normalizing p(beta)) " ){
        THEN( "returned value has as most a 1e-6 percent error" ){
            std::vector<double> p = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
            double tau = 0.5, delta = 1.0;
            equal( fsum( 1, p, tau, delta ), 39.387006 );

            p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
            delta = 2.0; tau = 0.5;
            equal( fsum( 1, p, tau, delta ), 29610795.321201 );

        } // THEN
    } // WHEN
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
        THEN( "returned value has as most a 1e-6 percent error" ){
            std::vector<double> p = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
            double tau = 0.5, delta = 0.1;
            equal( fsum( 0, p, tau, delta ), 0.035514341 );

            p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
            delta = 2.0; tau = 0.5;
            equal( fsum( 0, p, tau, delta ), 1444532.8400 );

        } // THEN
    } // WHEN
    WHEN( "n = 2 (used for effective temperature calculation) " ){
        THEN( "returned value has as most a 1e-6 percent error" ){
            std::vector<double> p = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
            double tau = 0.5, delta = 0.7;
            equal( fsum( 2, p, tau, delta ), 3.21945524 );

            p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
            delta = 2.0; tau = 0.5;
            equal( fsum( 2, p, tau, delta ), 612298146.17046 );

        } // THEN
    } // WHEN
} // TEST CASE

TEST_CASE( "normalize" ){
    GIVEN( "some vector p, spacing delta, normalizing value tbeta" ){
        std::vector<double> p = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
        double tbeta = 0.8, delta = 0.5;
        normalize( p, delta, tbeta );
        std::vector<double> correct= {2.62155807E-2, 5.24311614E-2, 
              7.86467421E-2, 0.104862322, 0.1310779, 0.15729348};
        THEN( "each value in vector has as most a 1e-6 percent error" ){
            equal_vec( p, correct );
        } // THEN

        p = {7.802837847e-3, 1.560567569e-2, 2.340851354e-2, 3.121135138e-2, 
             3.901418923e-2, 4.681702708e-2, 5.461986492e-2, 6.242270277e-2, 
             7.022554062e-2, 7.802837847e-2, 8.583121631e-2, 9.363405416e-2};
        delta = 0.2; tbeta = 2.0;
        normalize( p, delta, tbeta );
        correct = {5.2978935E-2, 0.10595787, 0.15893680, 0.21191574, 
                   0.26489467, 0.31787361, 0.37085254, 0.42383148, 
                   0.47681042, 0.52978935, 0.58276829, 0.63574722};
        THEN( "each value in vector has as most a 1e-6 percent error" ){
          equal_vec( p, correct );
        } //  THEN

    } // GIVEN 
    GIVEN( "two vectors where one is the first scaled by some constant" ){
        std::vector<double> p1 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
        std::vector<double> p2 = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
        double tbeta = 5.0, delta = 0.8;
        THEN( "outputs should be same, and both within 1e-6 to true value" ){
            normalize( p1, delta, tbeta );
            normalize( p2, delta, tbeta );
            std::vector<double> correct= {3.095826174E-2, 6.191652348E-2, 
              9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};
            equal_vec( p1, correct );
            equal_vec( p1, correct );
        } // THEN
    } // GIVEN
} // TEST CASE


