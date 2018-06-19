#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/310.h"

void checkVecHasZeros(std::vector<double>& v, int index_0, int index_1 ){
  for ( int i = index_0; i < index_1; ++i ){
    REQUIRE( 0.0 == Approx(v[i]).epsilon(1e-6) );
  }
}
void checkVecHasLin(std::vector<double>& v, int index_0, int index_1 ){
  for ( int i = index_0; i < index_1; ++i ){
    REQUIRE( (i+1) == Approx(v[i]).epsilon(1e-6) );
  }
}

void initializeLin(std::vector<double>& v){
  for ( size_t i = 0; i < v.size(); ++i ){ 
    v[i] = i + 1;
  }
}


TEST_CASE( "310" ){
  std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.97e-4, 3.5e-4, 4.2e-4, 5.06e-4,
    6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.03e-3, 
    2.277e-3, 2.6e-3, 3e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5e-3, 5.6e-3, 6.325e-3, 
    7.2e-3, 8.1e-3, 9.108e-3, 1e-2, 1.063e-2, 1.15e-2, 1.2397e-2, 1.33e-2, 
    1.417e-2, 1.5e-2, 1.6192e-2, 1.82e-2, 1.99e-2, 2.0493e-2, 2.15e-2, 2.28e-2,
    2.53e-2, 2.8e-2, 3.0613e-2, 3.38e-2, 3.65e-2, 3.95e-2, 4.2757e-2, 4.65e-2,
    5e-2, 5.6925e-2, 6.25e-2, 6.9e-2, 7.5e-2, 8.1972e-2, 9e-2, 9.6e-2, 0.1035,
    0.111573, 0.12, 0.128, 0.1355, 0.145728, 0.16, 0.172, 0.184437, 0.2, 0.2277,
    0.2510392, 0.2705304, 0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 
    0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.7, 0.78, 0.86, 0.95, 1.05, 
    1.16, 1.28, 1.42, 1.55, 1.70, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 
    3.42, 3.75, 4.07, 4.46, 4.90, 5.35, 5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 
    9.85, 10.0 };
    int ie = 0, jbeta = 0, iskip = 0, j = 0, lasym = 0, nbeta = 80;
    double enow = 0.0, temp = 296.0, bk = 8.617385e-5, break_val = 3000.0,
      therm = 2.53e-2, ep = 0.0, tev = 2.55074596e-2, az = 11.9, tevz = 2.53e-2,
      iinc = 2, lat = 1, az2 = 0.0, teff2 = 0.0, cliq = 0.0, 
      sb = 5.5348570016241778, sb2 = 0.0, teff = 6.1475562851499993E-2;

    std::vector<double> esi(95,0.0), xsi(95,0.0), ubar(118,0.0), x(20,0.0), 
      p2(118,0.0), p3(118,0.0);





  GIVEN( "temperature is lower than the break value of 3000" ){
    WHEN( "lasym is zero, vector entries are zero" ){
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, 
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "ie's incremented, jbeta = -#beta, vectors set to 0 (except esi)" ){
        REQUIRE( 1 == ie ); REQUIRE(   0 == iskip); 
        REQUIRE( 0 == j  ); REQUIRE( -80 == jbeta);

        REQUIRE( 1.0e-5 == Approx(  enow  ).epsilon(1e-6) );
        REQUIRE( 0.0    == Approx(  ep    ).epsilon(1e-6) );
        REQUIRE( 1.0e-5 == Approx( esi[0] ).epsilon(1e-6) );
  
        checkVecHasZeros( esi, 1, esi.size() );
        checkVecHasZeros( xsi, 0, xsi.size() );
        checkVecHasZeros( ubar, 0, ubar.size() );
      } // THEN
    } // WHEN


    WHEN( "lasym is some positive value" ){
      lasym = 1;
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, 
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "jbeta is set to equal 1" ){
        REQUIRE(1 == jbeta);
      } // THEN
    } // WHEN

    WHEN( "vectors contain nonzero entries, ie is nonzero index" ){
      ie = 3, jbeta = 10, iskip = 2, j = 8; enow = 0.5, ep = 1.0;

      initializeLin(esi); initializeLin(xsi); initializeLin(ubar); 
      initializeLin(x);   initializeLin(p2);  initializeLin(p3);
  
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2,
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "ie -> ie+1, and the ie+1 value in each vec is changed (except x)" ){
  
        REQUIRE( 4 == ie );     REQUIRE( -80 == jbeta );
        REQUIRE( 0  == iskip ); REQUIRE(   0 == j );
  
        REQUIRE( 3.5e-5 == Approx(enow).epsilon(1e-6) );
        REQUIRE( 0.0 == Approx(ep).epsilon(1e-6) );
  
        checkVecHasLin(esi,  0, 3 );
        checkVecHasLin(xsi,  0, 3 );
        checkVecHasLin(ubar, 0, 3 );
        REQUIRE( 3.5e-5 == Approx(esi[3]).epsilon(1e-6) );
        REQUIRE( 0.0    == Approx(xsi[3]).epsilon(1e-6) );
        REQUIRE( 0.0    == Approx(ubar[3]).epsilon(1e-6) );
        checkVecHasLin(esi,  4, esi.size() );
        checkVecHasLin(xsi,  4, xsi.size() );
        checkVecHasLin(ubar, 4, ubar.size() );
  
        REQUIRE( 0.0 == Approx(x[0]).epsilon(1e-6) );
        checkVecHasLin(x, 1, x.size() );
  
  
      } // THEN
    } // WHEN
  } // GIVEN




  GIVEN( "temperature is higher than the break value of 3000" ){
    temp = 3001;
    WHEN( "lasym is zero, vector entries are zero" ){
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, 
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "ie is incremented, jbeta = -#beta, vectors set to 0 (except esi)" ){
        REQUIRE( 1 == ie ); REQUIRE(   0 == iskip); 
        REQUIRE( 0 == j  ); REQUIRE( -80 == jbeta);

        REQUIRE( 1.0e-5 == Approx(  enow  ).epsilon(1e-6) );
        REQUIRE( 0.0    == Approx(  ep    ).epsilon(1e-6) );
        REQUIRE( 1.0e-5 == Approx( esi[0] ).epsilon(1e-6) );
  
        checkVecHasZeros( esi, 1, esi.size() );
        checkVecHasZeros( xsi, 0, xsi.size() );
        checkVecHasZeros( ubar, 0, ubar.size() );
      } // THEN
    } // WHEN


    WHEN( "lasym is some positive value" ){
      lasym = 1;
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, 
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "jbeta is set to equal 1" ){
        REQUIRE(1 == jbeta);
      } // THEN
    } // WHEN

    WHEN( "vectors contain nonzero entries, ie is nonzero index" ){
      ie = 3, jbeta = 10, iskip = 2, j = 8; enow = 0.5, ep = 1.0;

      initializeLin(esi); initializeLin(xsi); initializeLin(ubar); 
      initializeLin(x);   initializeLin(p2);  initializeLin(p3);
  
      do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2,
        p3, ep, jbeta, iskip, j, nbeta, lasym, x );

      THEN( "ie -> ie+1, and the ie+1 value in each vec is changed (except x)" ){
  
        REQUIRE( 4 == ie );     REQUIRE( -80 == jbeta );
        REQUIRE( 0  == iskip ); REQUIRE(   0 == j );
  
        REQUIRE( 4.32125e-5 == Approx(enow).epsilon(1e-6) );
        REQUIRE( 0.0 == Approx(ep).epsilon(1e-6) );
  
        checkVecHasLin(esi,  0, 3 );
        checkVecHasLin(xsi,  0, 3 );
        checkVecHasLin(ubar, 0, 3 );
        REQUIRE( 4.32125e-5 == Approx(esi[3]).epsilon(1e-6)  );
        REQUIRE( 0.0        == Approx(xsi[3]).epsilon(1e-6)  );
        REQUIRE( 0.0        == Approx(ubar[3]).epsilon(1e-6) );
        checkVecHasLin(esi,  4, esi.size() );
        checkVecHasLin(xsi,  4, xsi.size() );
        checkVecHasLin(ubar, 4, ubar.size() );
  
        REQUIRE( 0.0 == Approx(x[0]).epsilon(1e-6) );
        checkVecHasLin(x, 1, x.size() );
  
  
      } // THEN
    } // WHEN
  } // GIVEN

} // TEST CASE




