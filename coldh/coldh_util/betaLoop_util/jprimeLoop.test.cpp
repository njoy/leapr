#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../../catch.hpp"
#include "jprimeLoop.h"


void equal( double a, double b ){
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
      for ( auto i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }
      auto out1 = jPrime( total1, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, false );


      double total2 = 0.02; x = 0.85; swe = 0.87; pj = 0.48; y = 0.35;
      sex = {5,4,3,2,1,1,1,2,2,1,5};
      auto out2 = jPrime( total2, j, be, x, swe, pj, jj, bex, rdbex, sex, betan,
                         al, wt, tbart, y, nbx, false, false );

      THEN( "output matches expectation" ){
        equal( out1,   0.2353529421 );
        equal( total1, 4.2428657E-2 );
        equal( out2,   8.87488195E-2 );
        equal( total2, 8.65562642E-2 );
      } // THEN
    } // WHEN

    WHEN( "molecular translations are assumed to be free (free = true)" ){

    } // WHEN
  } // GIVEN


  GIVEN( "jprime loop for odd values" ){
    WHEN( "inputs" ){
    double total = 0.02, wt = 2.3, be = -1.2, al = 0.1, x = 0.85, 
           swo = 0.87, pj = 0.48, y = 0.35, tbart = 950.0;
    int j = 1, jj = 0, nbx = 10;
    std::vector<double> bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0 };
    std::vector<double> rdbex { 1.6666666, 3.3333333, 6.6666666, 20.0, 5.0, 
      20.0, 6.66666, 3.333333, 1.666666, 0.0, 0.0 };
    std::vector<double> sex { 5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.721415952, 
      2.2224546, 2.1952465, 1.5059710, 5.0 };
    std::vector<double> betan { 0.1, 0.15, 0.3, 0.6, 1.2 };

    double snlk = jPrime( total, j, be, x, swo, pj, jj, bex, 
        rdbex, sex, betan, al, wt, tbart, y, nbx, true, false );
    equal( snlk,  8.01757644 );
    equal( total, 1.62384370 );
    
    } // GIVEN

  WHEN( "other inputs" ){
    double total = 0.0, wt = 2.3, be = -1.2, al = 0.1, x = 0.8, 
           swo = 0.9, pj = 0.4, y = 0.3, tbart = 950.0;
    int j = 1, jj = 0, nbx = 10;
    std::vector<double> bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0 };
    std::vector<double> rdbex { 1.6666666, 3.3333333, 6.6666666, 20.0, 5.0, 
      20.0, 6.66666, 3.333333, 1.666666, 0.0, 0.0 };
    std::vector<double> betan { 0.1, 0.15, 0.3, 0.6, 1.2 };
    std::vector<double> sex(11,0.0);
    for ( auto i = 0; i < sex.size(); ++i ){
      sex[i] = i + 1;
    }
    double snlk = jPrime( total, j, be, x, swo, pj, jj, bex, 
        rdbex, sex, betan, al, wt, tbart, y, nbx, true, false );
    equal( snlk,  1.397417419 );
    equal( total, 1.397570948 );

  } // WHEN
  } // GIVEN

} // TEST CASE
