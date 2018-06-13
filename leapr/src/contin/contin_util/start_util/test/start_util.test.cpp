#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/start_util.h"


TEST_CASE( "fsum" ){
  std::vector<double> p;
  double tau, delta;

  WHEN( "n = 1 (used for normalizing p(beta)) " ){
    THEN( "returned value has as most a 1e-6 percent error" ){
      p = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
      tau = 0.5, delta = 1.0;
      REQUIRE( fsum(1,p,tau,delta) == Approx(39.387006).epsilon(1e-6) );

      p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
      delta = 2.0; tau = 0.5;
      REQUIRE( fsum(1,p,tau,delta) == Approx(29610795.32).epsilon(1e-6));

    } // THEN
  } // WHEN
  WHEN( "n = 0 (used for debye-waller coefficient) " ){
    THEN( "returned value has as most a 1e-6 percent error" ){
      p = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
      tau = 0.5, delta = 0.1;
      REQUIRE( fsum(0,p,tau,delta) == Approx(0.035514341).epsilon(1e-6) );

      p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
      delta = 2.0; tau = 0.5;
      REQUIRE( fsum(0,p,tau,delta) == Approx(1444532.8400).epsilon(1e-6) );

    } // THEN
  } // WHEN
  WHEN( "n = 2 (used for effective temperature calculation) " ){
    THEN( "returned value has as most a 1e-6 percent error" ){
      p = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
      tau = 0.5, delta = 0.7;
      REQUIRE( fsum(2,p,tau,delta) == Approx(3.21945524).epsilon(1e-6) );

      p = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
      delta = 2.0; tau = 0.5;
      REQUIRE( fsum(2,p,tau,delta) == Approx(612298146.17046).epsilon(1e-6) );

    } // THEN
  } // WHEN
} // TEST CASE

TEST_CASE( "normalize" ){
  std::vector<double> p;
  double tbeta, delta;
  GIVEN( "some vector p, spacing delta, normalizing value tbeta" ){

    p = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
    tbeta = 0.8, delta = 0.5;
    normalize( p, delta, tbeta );
    std::vector<double> correct= {2.62155807E-2, 5.24311614E-2, 
      7.86467421E-2, 0.104862322, 0.1310779, 0.15729348};

    THEN( "each value in vector has as most a 1e-6 percent error" ){
      REQUIRE( p.size() == correct.size() );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
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
      REQUIRE( p.size() == correct.size() );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
    } // THEN
  } // GIVEN 

  GIVEN( "two vectors where one is the first scaled by some constant" ){
    std::vector<double> p1 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12}, 
                        p2 = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    tbeta = 5.0, delta = 0.8;

    normalize( p1, delta, tbeta );
    normalize( p2, delta, tbeta );

    std::vector<double> correct= {3.095826174E-2, 6.191652348E-2, 
      9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};

    THEN( "outputs should be same, and both within 1e-6 to true value" ){
      REQUIRE( p1.size() == correct.size() );
      REQUIRE( p2.size() == correct.size() );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p1[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
        REQUIRE( p2[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
    } // THEN
  } // GIVEN
} // TEST CASE


