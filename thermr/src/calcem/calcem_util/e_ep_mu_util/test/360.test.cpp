#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/360.h"


void vecHasZeros( const std::vector<double>& v, size_t i1, size_t i2 ){
  for ( size_t i = i1; i < i2; ++i ){
    REQUIRE( 0.0 == Approx( v[i] ).epsilon(1e-6) );
  }
}

TEST_CASE( "360" ){
  int i = 4, j = 4, jmax = 55550, ie = 1, nll = 32673, nl = 9;
  std::vector<double> xsi(95,0.0), x(20,0.0), p2(118,0.0), p3(118,0.0), 
    p(4,0.0), ubar(118,0.0);
  x[0] = 1e-5; x[1] = 5e-5; x[2] = 2.5e-5;
  double ulast = 1.1, u2last = 1.2, u3last = 1.3, xlast = 2.3, ylast = 3.4,
         tolmin = 5.0e-6; 
  std::vector<std::vector<double>> y(65,std::vector<double>(20));
  for ( size_t i = 0; i < y.size(); ++i ){
    for ( size_t j = 0; j < y[0].size(); ++j ){
      y[i][j] = 0.01*(i+j);
    }
  }

  GIVEN( " " ){
    WHEN( "ie and i invoke the first entries of vectors and y, respectively" ){
      for ( size_t i = 0; i < y.size(); ++i ){
        for ( size_t j = 0; j < y[0].size(); ++j ){
          y[i][j] = 0.01*(i+j);
        }
      }

      AND_WHEN( "j == 3 and xsi(ie) < tolmin" ){
        THEN( "j is incremented by 1" ){
          do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, u3last, tolmin, 
            y, p2, p3, nll, nl, p, ubar, ie, i );
          REQUIRE( 5 == j );

          REQUIRE( -3.9445     == Approx(xsi[0]).epsilon(1e-6)  );
          REQUIRE( -1.2675875  == Approx(ubar[0]).epsilon(1e-6) );
          REQUIRE( -1.36306828 == Approx(p2[0]).epsilon(1e-6)   );
          REQUIRE( -1.49116523 == Approx(p3[0]).epsilon(1e-6)   );
          vecHasZeros( xsi,  1, xsi.size()  );
          vecHasZeros( ubar, 1, ubar.size() );
          vecHasZeros( p2,   1, p2.size()   );
          vecHasZeros( p3,   1, p3.size()   );
 
        
          for ( size_t i = 1; i < xsi.size(); ++i ){ 
            REQUIRE( 0.0 == Approx(xsi[i]).epsilon(1e-6) ); 
          }

        } // THEN
      } // AND WHEN
      AND_WHEN( "j == 3 and xsi(ie) < tolmin" ){
        j = 2;
        THEN( "j is set to be 2" ){
          do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, u3last, tolmin, 
            y, p2, p3, nll, nl, p, ubar, ie, i );
          REQUIRE( j == 2 );
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "ie and i invoke later entries of vectors and y, respectively" ){
      for ( size_t i = 0; i < y.size(); ++i ){
        for ( size_t j = 0; j < y[0].size(); ++j ){
          y[i][j] = 0.2*(i+j);
        }
      }
      xlast = 0.3; ylast = 0.4; ulast = 0.1; u2last = 0.2; u3last = 0.3;
      ie = 3; i = 4;
      do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, u3last, tolmin, y, 
        p2, p3, nll, nl, p, ubar, ie, i );
      REQUIRE( -0.15 == Approx(xsi[2] ).epsilon(1e-6) );
      REQUIRE( -0.15 == Approx(ubar[2]).epsilon(1e-6) );
      REQUIRE( -0.3171 == Approx(p2[2]  ).epsilon(1e-6) );
      REQUIRE( -0.8145 == Approx(p3[2]  ).epsilon(1e-6) );

      vecHasZeros( xsi,  0, 2 ); vecHasZeros( xsi,  3, xsi.size()  );
      vecHasZeros( ubar, 0, 2 ); vecHasZeros( ubar, 3, ubar.size() ); 
      vecHasZeros( p2,   0, 2 ); vecHasZeros( p2,   3, p2.size()   ); 
      vecHasZeros( p3,   0, 2 ); vecHasZeros( p3,   3, p3.size()   ); 
 
    } // WHEN
  } // GIVEN



  GIVEN( "iinvalid (out of range) value of j" ){
    j = 55551;
    THEN( "an exception is thrown" ){
      REQUIRE_THROWS( do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, 
        u3last, tolmin, y, p2, p3, nll, nl, p, ubar, ie, i ) );
    } // THEN
  } // GIVEN
} // TEST CASE
