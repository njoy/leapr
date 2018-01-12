#include <iostream>
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop.h"


void equal1( double a, double b ){
  if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}


TEST_CASE( "jprime loop" ){
  GIVEN( "jprime loop is over even values" ){

    std::vector<double> 
    bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 },
    rdbex { 1.6667, 3.33333, 6.6667, 20, 5, 20, 6.6667, 3.333333, 1.6667, 0, 0 },
    betan { 0.1, 0.15, 0.3, 0.6, 1.2 },
    sex(11,0.0);

    int j = 1, ifree = 0, jj = 0, nbx = 10;

    double wt = 2.3, be = -1.2, al = 0.1, x = 0.8, swe = 0.9, 
             pj = 0.4, y = 0.3, tbart = 950.0;

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
    
      double total1 = 0.0;
      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }
      auto out1 = jPrime( total1, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, false );


      double total2 = 0.02; x = 0.85; swe = 0.87; pj = 0.48; y = 0.35;
      sex = {5,4,3,2,1,1,1,2,2,1,5};
      auto out2 = jPrime( total2, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, false );

      THEN( "output matches expectation" ){
        equal1( out1,   0.2353529421 );
        equal1( total1, 4.2428657E-2 );
        equal1( out2,   8.87488195E-2 );
        equal1( total2, 8.65562642E-2 );
      } // THEN
    } // WHEN

    WHEN( "molecular translations are assumed to be free (free = true)" ){
      double total1 = 3; j = 4; be = 1.23; x = 0.1; swe = 0.2; pj = 2.4; 
             jj = 2; al = 2.0; wt = 0.05; y = 0.3;

      auto out1 = jPrime( total1, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, true );

      double total2 = 5.0; jj = 0;

      auto out2 = jPrime( total2, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, true );

      THEN( "output matches expectation" ){
        equal1( total1, 3.0 );
        equal1( out1,   1.9988998E-2 );
        equal1( total2, 6.8634280542 );
        equal1( out2,   1.9988998E-2 );

      } // THEN

    } // WHEN
  } // GIVEN


  GIVEN( "jprime loop for odd values" ){

    std::vector<double> 
    bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 },
    rdbex{ 1.666667, 3.33333, 6.6667, 20, 5, 20, 6.6667, 3.3333, 1.66667, 0, 0 },
    betan { 0.1, 0.15, 0.3, 0.6, 1.2 },
    sex(11,0.0);

    double wt = 2.3, be = -1.2, al = 0.1, x = 0.85, 
           swo = 0.87, pj = 0.48, y = 0.35, tbart = 950.0;

    int j = 1, jj = 0, nbx = 10;

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      sex = { 5, 4, 3, 2, 1, 1, 1.72141595, 2.222455, 2.195246, 1.505971, 5 };
      double total1 = 0.02;
      double out1 = jPrime( total1, j, be, x, swo, pj, jj, bex, 
          rdbex, sex, betan, al, wt, tbart, y, nbx, true, false );

      THEN( "output matches expectation" ){
        equal1( out1,   8.01757644 );
        equal1( total1, 1.62384370 );
      } // THEN

      double total2 = 0.0; x = 0.8; swo = 0.9; pj = 0.4; y = 0.3;
      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }
      double out2= jPrime( total2, j, be, x, swo, pj, jj, bex, 
        rdbex, sex, betan, al, wt, tbart, y, nbx, true, false );

      THEN( "output matches expectation" ){
        equal1( out2,   1.397417419 );
        equal1( total2, 1.397570948 );
      } // THEN


    } // WHEN

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = 2 * (i + 1); }
      double total = 0.02; wt = 2.05; y = 0.35; swo = 0.47; total = 5.0;
      double out = jPrime( total, j, be, x, swo, pj, jj, bex, 
          rdbex, sex, betan, al, wt, tbart, y, nbx, true, true );

      THEN( "output matches expectation" ){
        equal1( out,   0.16137096895 );
        equal1( total, 5.86644429347 );
      } // THEN
    } // WHEN

  } // GIVEN

} // TEST CASE
