#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sint.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}



TEST_CASE( "sint" ){
  GIVEN( "inputs" ){
    double x = -0.1, alpha = 0.1, wt = 2.3, tbart = 407.4545311;
    int b = 4, nbx = 10;
    std::vector<double> bex {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0};
    std::vector<double> rdbex  {1.66666667, 3.33333333, 6.666666667, 20.0, 5.0,
      20.0, 6.66666667, 3.33333333, 1.66666667, 0.0, 0.0};
    std::vector<double> betan {0.1, 0.15, 0.3, 0.6, 1.2};
    std::vector<double> sex  {5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.7214159, 
      2.2224546, 2.1952465, 1.5059710, 0.0};
  

    auto sintOut = sint( x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx );
    equal( sintOut, 1.0 );


    WHEN( "OPTION A" ){
      THEN( "correct values are output" ){
        x = -1.201; alpha = -0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 0.0 );
        x = 1.201; alpha = -1e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 0.0 );
      } // THEN
    } // WHEN

    WHEN( "OPTION B" ){
      THEN( "correct values are output" ){
        x = 1.201; alpha = 0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 2.548608E-4 );

        x = -1.201; alpha = 1e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 1.654202E-16 );
      
      } // THEN
    } // WHEN

    WHEN( "OPTION C" ){
      THEN( "correct values are output" ){
        x = -0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 1.0 );
        x = -0.15;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 2.0 );

        x = -0.3;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 3.0 );


      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
