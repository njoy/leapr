#define CATCH_CONFIG_MAIN
#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/hexLatticeFactors_util/hexLatticeFactorsHelper.h"


TEST_CASE( "Function to Compute Hexagonal Lattice Factors" ){
  int k;
  double tsq, tsqx, f;
  std::vector<double> b ( 60000, 0.0 );

  k = 0; 

  GIVEN( "k is not positive or tsq <= tsqx" ){
    tsq = 0.1; tsqx = 9.6;
    THEN( "b's first two entries are populted with tsq and f, respectively" ){
 
      std::vector<double> fVec {3.0, 100, 100}, tsqVec {0.1, 0.1, 50};
      std::vector<std::tuple<double,double,double>> 
        good {{0.1,3.0,1},{0.1,100,1},{50,100,1}};

      for ( size_t j = 0; j < 3; ++j ){
        tsq = tsqVec[j]; f = fVec[j]; k = 0;
        hexLatticeFactorsHelper( k, tsq, tsqx, b, f );
        REQUIRE( b[0] == Approx(std::get<0>(good[j])).epsilon(1e-6) );
        REQUIRE( b[1] == Approx(std::get<1>(good[j])).epsilon(1e-6) );
        REQUIRE( k    == Approx(std::get<2>(good[j])).epsilon(1e-6) );
        for ( size_t i = 2; i < b.size(); ++i ){ REQUIRE( b[i] == 0 ); }
      } 
    } // THEN
  } // GIVEN

  GIVEN( "k is positive and tsq > tsqx" ){
    WHEN( "loop is able to finish" ){
      THEN( "we're in the 2st situation, f goes in b (only 1 val added)" ){

        std::vector<double> fVec {4.5e-2, 1.5, 4.5e-1}, 
		tsqVec {30.1,10.1,60.1}, tsqxVec {25, 5, 15 }, 
		kVec {2, 1, 2};
	std::vector<std::vector<double>> bAnswer { {10,20,30,40.045,0,0}, 
          {10,21.5,30,40,0,0}, {10,20,30,40,60.1,0.45} };

        for ( size_t j = 0; j < 3; ++j ){
          tsq = tsqVec[j]; tsqx = tsqxVec[j]; f = fVec[j]; k = kVec[j];
          b[0] = 10; b[1] = 20; b[2] = 30; b[3] = 40;
          hexLatticeFactorsHelper( k, tsq, tsqx, b, f );

          for ( size_t i = 0; i < 6; ++i ){ 
            REQUIRE( b[i] == Approx(bAnswer[j][i]).epsilon(1e-6) );
          }
          for ( size_t i = 6; i < b.size(); ++i ){ REQUIRE( b[i] == 0 ); }
        } 
      } // THEN
    } // WHEN
    
    WHEN( "loop isn't able to finish" ){
      THEN( "we're in the second situation" ){
        tsq = 700, tsqx = 500;
        k = 1; f = 4.5e-2;
        hexLatticeFactorsHelper( k, tsq, tsqx, b, f );

        REQUIRE( b[0] == Approx(0).epsilon(1e-6) ); 
	REQUIRE( b[1] == Approx(0).epsilon(1e-6) ); 
	REQUIRE( b[2] == Approx(700).epsilon(1e-6) ); 
        REQUIRE( b[3] == Approx(4.5e-2).epsilon(1e-6) ); 
	REQUIRE( k == Approx(2).epsilon(1e-6) ); 
        for ( size_t i = 4; i < b.size(); ++i ){ REQUIRE( b[i] == 0 ); }

        tsq = 63, tsqx = 50;
        k = 1; f = 1.5;
        b[0] = 10; b[1] = 20; b[2] = 30; b[3] = 40; b[4] = 50; b[5] = 60;
        hexLatticeFactorsHelper( k, tsq, tsqx, b, f );

        REQUIRE( b[0] == Approx(10).epsilon(1e-6) ); 
	REQUIRE( b[1] == Approx(20).epsilon(1e-6) ); 
	REQUIRE( b[2] == Approx(63).epsilon(1e-6) ); 
        REQUIRE( b[3] == Approx(1.5).epsilon(1e-6) ); 
	REQUIRE( b[4] == Approx(50).epsilon(1e-6) ); 
	REQUIRE( b[5] == Approx(60).epsilon(1e-6) );
	REQUIRE( k == Approx(2).epsilon(1e-6) ); 
        for ( size_t i = 6; i < b.size(); ++i ){ REQUIRE( b[i] == 0 ); }


      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE

