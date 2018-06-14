#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coh/coh_util/upstk.h"

void checkStk( const std::vector<std::vector<double>>& stk, 
  const std::vector<double>& stkVals ){
  // stkVals are only the nonzero values in the desired output. In the 
  // for loop down there, we want all values beyond stkVals's range to be 0

  for ( int i = 0; i < stk[0].size(); ++i ){
    for ( int j = 0; j < stk.size(); ++j ){
      if ( stk.size()*i < stkVals.size() ){
        REQUIRE( stk[j][i] == Approx(stkVals[stk.size()*i+j]).epsilon(1e-6) );
      } 
      else {
        REQUIRE( stk[j][i] == Approx(0.0).epsilon(1e-6) );
      }
    }
  }
} 

auto initializeStk( int nx, int maxIter ){
  std::vector<std::vector<double>> stk ( nx, std::vector<double> (20,0.0) );
  for ( int i = 0; i < maxIter; ++i ){
    for ( int j = 0; j < nx; ++j ){
      stk[j][i] = i*0.1*nx + 0.1*(j + 1);
    }
  }
  return stk;
}


TEST_CASE( "uptsk" ){
  int nl, nx, i = 1;
  std::vector<double> s(6), stkVals(18);
  double e;

  GIVEN( "legendre order nl of 1" ){
    nl = 1;
    e = 4.0e-4;

    WHEN( "cycle length is small (nx=2)" ){
      nx = 2;
      s = {0,0,0,0,0,0};

      AND_WHEN( "stk is initialized to nearly all zeros" ){

        std::vector<std::vector<double>> stk (nx,std::vector<double> (20,0.0));
        stk[0][0] = 4.5e-4;

        upstk( e, s, stk, nl, nx, i );

        THEN( "output is correct" ){
          stkVals = { 4.0e-4, 0, 4.5e-4, 0 };
          checkStk( stk, stkVals );
        } // THEN

      } // AND WHEN
  
      AND_WHEN( "stk is initialized to some nonzero values" ){

        auto stk = initializeStk( nx, 5 );

        upstk( e, s, stk, nl, nx, i );

        stkVals = {4e-4, 0, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

        THEN( "output is correct" ){
          checkStk( stk, stkVals );
        } // THEN

      } // AND WHEN

    } // WHEN

    WHEN( "cycle length is moderate (nx=4)" ){
      s[0] = 1.5e3;
      nx = 4;
      e = 1.8e-3;

      auto stk = initializeStk( nx, 3 );

      upstk( e, s, stk, nl, nx, i );

      stkVals = { 1.8e-3, 1500, 0.3, 0.4, 0.1, 0.2, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2 };

      THEN( "output is correct" ){
        checkStk( stk, stkVals );
      } // THEN

    } // WHEN

  } // GIVEN

  GIVEN( "legendre order nl of 3" ){
    s = { 1.5e3, 2.5e3, 3.5e3, 4.5e3, 5.5e3, 6.5e3 };
    e = 1.8e-3;
    nl = 3;

    WHEN( "cycle length is large (nx=6)" ){
      nx = 6;

      auto stk = initializeStk( nx, 3 );

      upstk( e, s, stk, nl, nx, i );

      stkVals = { 1.8e-3, 1500, 2500, 3500, 0.5, 0.6, 0.1, 0.2, 0.3, 0.4, 1.1, 
        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 };
      
      THEN( "output is correct" ){
        checkStk( stk, stkVals );
      } // THEN

    } // WHEN

    WHEN( "cycle length is moderate (nx=4)" ){
      nx = 4;

      auto stk = initializeStk( nx, 3 );

      upstk( e, s, stk, nl, nx, i );

      stkVals = { 1.8e-3, 1500, 2500, 3500, 0.1, 0.2, 0.3, 0.4, 0.9, 1.0, 1.1, 1.2 };

      THEN( "output is correct" ){
          checkStk( stk, stkVals );
      } // THEN

    } // WHEN

  } // GIVEN

} // TEST CASE



