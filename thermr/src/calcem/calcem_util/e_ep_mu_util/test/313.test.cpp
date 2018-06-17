#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/313.h"

TEST_CASE( "313" ){
  int jbeta = 10, lat = 1, iskip = 0;
  double ep = 2.5, enow = 1.3, tev = 2.5e-2, tevz = 2.53e-2;
  std::vector<double> beta (80), x(20);
  for ( size_t i = 0; i < 80; ++i ){ beta[i] = 0.20*i + 0.05; }

  GIVEN( "negative jbeta" ){
    jbeta = -13; lat = 1; enow = -0.95;
    WHEN( "lat = 1" ){
      AND_WHEN( "x vector entries are relatively small" ){
        jbeta = -3;
        for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.03*i - 5.0; }
        THEN( "single iteration is necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( -3 == jbeta );
          REQUIRE( -0.96138499 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
      AND_WHEN( "x vector entries are relatively large" ){
        for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.1*i - 0.9; }
        THEN( "multiple iterations are necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 31 == jbeta );
          REQUIRE( -0.79693498 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "lat != 1" ){
      lat = 0;
      AND_WHEN( "x vector entries are relatively small" ){
        jbeta = -3;
        for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.03*i - 5.0; }
        THEN( "single iteration is necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( -3 == jbeta );
          REQUIRE( -0.96124999 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
      AND_WHEN( "x vector entries are relatively large" ){
        THEN( "multiple iterations are necessary" ){
          for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.1*i - 0.9; }
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 31 == jbeta );
          REQUIRE( -0.79874998 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
          jbeta = -13; enow = -0.95; 

          for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.03*(i+1) - 1.0; }
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 3 == jbeta );
          REQUIRE( -0.93874999 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN

      } // AND WHEN
    } // WHEN
    WHEN( "E' equals E_now" ){
      lat = 1; jbeta = -3; enow = -0.95;
      for ( size_t i = 0; i < 80; ++i ){ beta[i] = 0.2*i - 0.2; }
      for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.1*i - 0.9; }
      THEN( "E' is rounded slightly up" ){
        do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 32 == jbeta );
          REQUIRE( -0.79819998 == Approx(ep).epsilon(1e-6) );
          REQUIRE( -0.94999998 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 1 == iskip );  

      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "positive or zero jbeta" ){
    WHEN( "lat = 1" ){
      AND_WHEN( "x vector entries are relatively small" ){
        for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.03*i + 0.03; }
        THEN( "single iteration is necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 10 == jbeta );
          REQUIRE( 1.346805 == Approx(ep).epsilon(1e-6) );
          REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
      AND_WHEN( "x vector entries are relatively large" ){
        for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.7*i + 0.7; }
        THEN( "multiple iterations are necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 21 == jbeta );
          REQUIRE( 1.402465 == Approx(ep).epsilon(1e-6) );
          REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "lat != 1" ){
      lat = 0;
      AND_WHEN( "x vector entries are relatively small" ){
        for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.03*i + 0.03; }
        THEN( "single iteration is necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 10 == jbeta );
          REQUIRE( 1.34625 == Approx(ep).epsilon(1e-6) );
          REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
      AND_WHEN( "x vector entries are relatively large" ){
        for ( size_t i = 0; i < 20; ++i ){ x[i] = 0.7*i + 0.7; }
        THEN( "multiple iterations are necessary" ){
          do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 21 == jbeta );
          REQUIRE( 1.40125 == Approx(ep).epsilon(1e-6) );
          REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 0 == iskip );  
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "E' equals E_now" ){
      lat = 0;
      for ( size_t i = 0; i < 80; ++i ){ beta[i] = 0.2*i - 5.0; }
      for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.7*i + 0.7; }
      THEN( "E' is rounded slightly up" ){
        do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
          REQUIRE( 46 == jbeta );
          REQUIRE( 1.4 == Approx(ep).epsilon(1e-6) );
          REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
          REQUIRE( 1 == iskip );  

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
