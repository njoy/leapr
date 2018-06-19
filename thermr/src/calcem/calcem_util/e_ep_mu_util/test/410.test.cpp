#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/410.h"

TEST_CASE( "410" ){
  GIVEN( " " ){
    int i = 2, nl = 4;
    double xm = 5.0e-5;
    std::vector<double> x (20,0.0), yt (65,0.0);
    std::vector<std::vector<double>> y (65,std::vector<double>(20)); 
    for ( size_t i = 0; i < y.size(); ++i ){
      for ( size_t j = 0; j < y[0].size(); ++j ){
        y[i][j] = 0.01*(i+j);
      }
    }
    yt[0] =  0.3; yt[1] = -0.9; yt[2] = -0.7; yt[3] = -0.5; yt[4] = -0.4; 
    yt[5] = -0.2; yt[6] = 0.08; yt[7] =  0.2; yt[8] =  0.6;

    x[0] = 1e-5;
    WHEN( "lat = 1" ){
      do410(i,x,xm,nl,y,yt);
      std::vector<std::tuple<double,double,double,double,double>> output { 
        { 0.00, 0.30, 0.01, 0.03, 0.04 }, { 0.01,-0.90, 0.02, 0.04, 0.05 },
        { 0.02,-0.70, 0.03, 0.05, 0.06 }, { 0.03,-0.50, 0.04, 0.06, 0.07 },
        { 0.04, 0.05, 0.06, 0.07, 0.08 } };
      THEN( "E' is rounded slightly up" ){
        for ( size_t i = 0; i < output.size(); ++i ){
          REQUIRE( std::get<0>(output[i]) == Approx(y[i][0]).epsilon(1e-6) );
          REQUIRE( std::get<1>(output[i]) == Approx(y[i][1]).epsilon(1e-6) );
          REQUIRE( std::get<2>(output[i]) == Approx(y[i][2]).epsilon(1e-6) );
          REQUIRE( std::get<3>(output[i]) == Approx(y[i][3]).epsilon(1e-6) );
          REQUIRE( std::get<4>(output[i]) == Approx(y[i][4]).epsilon(1e-6) );
        }


        std::vector<double> correctX { 1.0e-5, 5.0e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        for ( size_t i = 0; i < x.size(); ++i ){ 
          REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) ); 
        }

        REQUIRE( 3 == i );


      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
