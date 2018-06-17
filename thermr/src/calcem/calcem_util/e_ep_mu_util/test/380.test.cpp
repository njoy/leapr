#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/380.h"

TEST_CASE( "380" ){
  GIVEN( " " ){
    int i = 4, nl = 9, j = 66, nll, jnz = 0;
    double em9 = 1.0e-9, xlast, ylast, ulast, u2last, u3last;
    std::vector<double> p (4,0.0), x(20,0.0), scr(500000);
    std::vector<std::vector<double>> y (65,std::vector<double>(20)); 
    for ( int i = 0; i < y.size(); ++i ){
      for ( size_t j = 0; j < y[0].size(); ++j ){
        y[i][j] = 0.01*(i+j);
      }
    }
    for ( int i = 0; i < scr.size(); ++i ){
      scr[i] = 1.2*(i%3)+0.01*i;
    }
    x[0] = 1.0e-5; x[1] = 5.0e-6; x[2] = 2.5e-6;


    WHEN( "" ){
      THEN( "scr is not changed" ){
        do380( i, j, nl, nll, scr, x, y, em9, xlast, ylast, jnz, ulast, u2last, u3last, p );
        REQUIRE( 3 == i );
        REQUIRE( 0.0     == Approx(xlast).epsilon(1e-6) );
        REQUIRE( 3.0e-2  == Approx(ylast).epsilon(1e-6) );
        REQUIRE( 2.25e-3 == Approx(ulast).epsilon(1e-6) );
        REQUIRE( -1.4723249e-2 == Approx(u2last).epsilon(1e-6) );
        REQUIRE( -3.3344998e-3 == Approx(u3last).epsilon(1e-6) );
        std::vector<double> scr_645_to_670 {6.45, 7.66, 8.87, 6.48, 7.69, 8.9, 
          6.51, 7.72, 8.93, 6.54, 7.75, 0.0, 3.0e-2, 4.0e-2, 5.0e-2, 6.0e-2, 
          7.0e-2, 8.0e-2, 9.0e-2, 0.1, 0.11, 6.66, 7.87, 9.08, 6.69, 7.9 };
        for ( int i = 0; i < 645; ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }
        for ( int i = 0; i < scr_645_to_670.size(); ++i ){
          REQUIRE( scr_645_to_670[i] == Approx(scr[i+645]).epsilon(1e-6) );
        }
        for ( int i = 670; i < scr.size(); ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }


      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( " " ){
    int i = 4, nl = 9, j = 30, nll, jnz = 0;
    double em9 = 1.0e-9, xlast, ylast, ulast, u2last, u3last;
    std::vector<double> p (4,0.0), x(20,0.0), scr(500000);
    std::vector<std::vector<double>> y (65,std::vector<double>(20)); 
    for ( int i = 0; i < y.size(); ++i ){
      for ( size_t j = 0; j < y[0].size(); ++j ){
        y[i][j] = 0.2*(i+j);
      }
    }
    for ( int i = 0; i < scr.size(); ++i ){
      scr[i] = 1.2*(i%3)+0.01*i;
    }
    x[0] = 1.0e-5; x[1] = 5.0e-6; x[2] = 2.5e-6;


    WHEN( "scr value is greater than 1" ){
      THEN( "overflow is fixed, plus warning" ){
        do380( i, j, nl, nll, scr, x, y, em9, xlast, ylast, jnz, ulast, u2last, u3last, p );
        REQUIRE( 3 == i );
        REQUIRE( 0.0 == Approx(xlast).epsilon(1e-6) );
        REQUIRE( 0.6 == Approx(ylast).epsilon(1e-6) );
        REQUIRE( 0.9 == Approx(ulast).epsilon(1e-6) );
        REQUIRE( 1.914 == Approx(u2last).epsilon(1e-6) );
        REQUIRE( 5.13  == Approx(u3last).epsilon(1e-6) );

        std::vector<double> scr_295_to_306 { 4.15, 0.0, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.06 };
          
        for ( int i = 0; i < 295; ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }
        for ( int i = 0; i < scr_295_to_306.size(); ++i ){
          REQUIRE( scr_295_to_306[i] == Approx(scr[i+295]).epsilon(1e-6) );
        }
        for ( int i = 306; i < scr.size(); ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }


      } // THEN
    } // WHEN

    WHEN( "scr value is lesser than 1" ){
    for ( int i = 0; i < y.size(); ++i ){
      for ( size_t j = 0; j < y[0].size(); ++j ){
        y[i][j] = -0.2*(i+j);
      }
    }

      j = 3;
      THEN( "overflow is fixed, plus warning" ){
        do380( i, j, nl, nll, scr, x, y, em9, xlast, ylast, jnz, ulast, u2last, u3last, p );
        REQUIRE( 3 == i );
        REQUIRE( 0.0 == Approx(xlast).epsilon(1e-6) );
        REQUIRE( -0.6 == Approx(ylast).epsilon(1e-6) );
        REQUIRE( 0.9 == Approx(ulast).epsilon(1e-6) );
        REQUIRE( -1.914 == Approx(u2last).epsilon(1e-6) );
        REQUIRE( 5.13  == Approx(u3last).epsilon(1e-6) );

        std::vector<double> scr_25_to_39 { 1.45, 0.0, -0.6, -0.8, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.36, 1.57, 2.78, 0.39 };
          
        for ( int i = 0; i < 25; ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }
        for ( int i = 0; i < scr_25_to_39.size(); ++i ){
          REQUIRE( scr_25_to_39[i] == Approx(scr[i+25]).epsilon(1e-6) );
        }
        for ( int i = 39; i < scr.size(); ++i ){
          REQUIRE( 1.2*(i%3)+0.01*i == Approx(scr[i]).epsilon(1e-6) );
        }


      } // THEN
    } // WHEN

  } // GIVEN

} // TEST CASE
