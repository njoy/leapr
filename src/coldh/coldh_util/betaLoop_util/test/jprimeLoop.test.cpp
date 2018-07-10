#include <iostream>
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop.h"


TEST_CASE( "jprime loop" ){
  std::vector<double> bex(11), rdbex(11), betan(5), sex(11);
  int j, jj, nbx;
  double wt, be, al, x, swe, pj, y, tbart, total, out, swo;

  bex   = {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0};
  rdbex = {1.666667, 3.3333, 6.6667, 20, 5, 20, 6.6667, 3.33333, 1.66667, 0, 0};
  betan = {0.1, 0.15, 0.3, 0.6, 1.2};

  GIVEN( "jprime loop is over even values" ){
    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      wt = 2.3; be = -1.2; al = 0.1; tbart = 950;
      AND_WHEN( "starting total is 0" ){
        j = 1; jj = 0; nbx = 10;
        x = 0.8; swe = 0.9; pj = 0.4; y = 0.3; 
	total = 0;
        for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }

        out = jPrime( total, j, be, x, swe, pj, jj, bex, rdbex, sex, betan, al,
          wt, tbart, y, nbx, false, false );

        THEN( "output matches expectation" ){
          REQUIRE( out   == Approx(0.2353529421).epsilon(1e-6) );
          REQUIRE( total == Approx(4.2428657E-2).epsilon(1e-6) );
        } // THEN
      } // AND WHEN

      AND_WHEN( "starting total and vector are not 0" ){
        j = 1; jj = 0; nbx = 10;
        x = 0.85; swe = 0.87; pj = 0.48; y = 0.35;

        total = 0.02; 
        sex = {5,4,3,2,1,1,1,2,2,1,5};

        out = jPrime( total, j, be, x, swe, pj, jj, bex, rdbex, sex, betan, al,
          wt, tbart, y, nbx, false, false );

        THEN( "output matches expectation" ){
          REQUIRE( out   == Approx(8.87488195E-2).epsilon(1e-6) );
          REQUIRE( total == Approx(8.65562642E-2).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
    } // WHEN


    WHEN( "molecular translations are assumed to be free (free = true)" ){

      AND_WHEN( "starting total is 0" ){
        total = 3; j = 4; be = 1.23; x = 0.1; swe = 0.2; pj = 2.4; 
          jj = 2; al = 2.0; wt = 0.05; y = 0.3;

        out = jPrime( total, j, be, x, swe, pj, jj, bex, rdbex, sex, betan, al,
          wt, tbart, y, nbx, false, true );


        THEN( "output matches expectation" ){
          REQUIRE( total == Approx(3.0).epsilon(1e-6) );
          REQUIRE( out == Approx(1.9988998E-2).epsilon(1e-6) );
        } // THEN
      } // AND WHEN

      AND_WHEN( "starting total is not 0" ){
        total = 5; j = 4; be = 1.23; x = 0.1; swe = 0.2; pj = 2.4; 
          jj = 0; al = 2.0; wt = 0.05; y = 0.3;

        out = jPrime( total, j, be, x, swe, pj, jj, bex, rdbex, sex, betan, al,
          wt, tbart, y, nbx, false, true );

        THEN( "output matches expectation" ){
          REQUIRE( total == Approx(6.8634280542).epsilon(1e-6) );
          REQUIRE( out == Approx(1.9988998E-2).epsilon(1e-6) );

        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN


  GIVEN( "jprime loop for odd values" ){

    wt = 2.3; be = -1.2; al = 0.1; tbart = 950.0;
    j = 1, jj = 0, nbx = 10;

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      AND_WHEN( "total starts at zero" ){
      
        total = 0.0; x = 0.8; swo = 0.9; pj = 0.4; y = 0.3;

        for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }
        out= jPrime( total, j, be, x, swo, pj, jj, bex, 
        rdbex, sex, betan, al, wt, tbart, y, nbx, true, false );

        THEN( "output matches expectation" ){
          REQUIRE( out == Approx(1.397417419).epsilon(1e-6) );
          REQUIRE( total == Approx(1.397570948).epsilon(1e-6) );
        } // THEN
      } // AND WHEN

      AND_WHEN( "start at not zero" ){
        sex = { 5, 4, 3, 2, 1, 1, 1.72141595, 2.222455, 2.195246, 1.505971, 5 };
        total = 0.02;
	x = 0.85; swo = 0.87; y = 0.35; pj = 0.48;
        out = jPrime( total, j, be, x, swo, pj, jj, bex, rdbex, sex, betan, al, 
          wt, tbart, y, nbx, true, false );

        THEN( "output matches expectation" ){
          REQUIRE( out == Approx(8.01757644).epsilon(1e-6) );
          REQUIRE( total == Approx(1.62384370).epsilon(1e-6) );
        } // THEN
      } // AND WHEN
    } // WHEN

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = 2 * (i + 1); }
      total = 0.02; wt = 2.05; y = 0.35; swo = 0.47; total = 5.0;
      out = jPrime( total, j, be, x, swo, pj, jj, bex, rdbex, sex, betan, al, 
        wt, tbart, y, nbx, true, true );

      THEN( "output matches expectation" ){
        REQUIRE( out == Approx(0.16137096895).epsilon(1e-6) );
        REQUIRE( total == Approx(5.86644429347).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
