
#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/formf.h"


TEST_CASE( "formf" ){
  using std::get;
  GIVEN( "graphite input" ){
    REQUIRE( Approx(formf(1,0,0,0)).epsilon(1e-6) == 4 );
    REQUIRE( Approx(formf(1,1,0,0)).epsilon(1e-6) == 0.25 );
    REQUIRE( Approx(formf(1,0,1,0)).epsilon(1e-6) == 0.25 );
    REQUIRE( Approx(formf(1,0,0,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(1,1,1,0)).epsilon(1e-6) == 4 );
    REQUIRE( Approx(formf(1,1,0,1)).epsilon(1e-6) == 0.75 );
    REQUIRE( Approx(formf(1,0,1,1)).epsilon(1e-6) == 0.75 );
    REQUIRE( Approx(formf(1,1,1,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(1,1,2,3)).epsilon(1e-6) == 0.75 );
    REQUIRE( Approx(formf(1,1,0,3)).epsilon(1e-6) == 0.75 );
    REQUIRE( Approx(formf(1,1,4,0)).epsilon(1e-6) == 4 );
    REQUIRE( Approx(formf(1,0,4,8)).epsilon(1e-6) == 0.25 );
  } // GIVEN

  GIVEN( "beryllium input" ){
    REQUIRE( Approx(formf(2,0,0,0)).epsilon(1e-6) == 2 );
    REQUIRE( Approx(formf(2,1,0,0)).epsilon(1e-6) == 0.5 );
    REQUIRE( Approx(formf(2,0,1,0)).epsilon(1e-6) == 0.5 );
    REQUIRE( Approx(formf(2,0,0,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(2,1,1,0)).epsilon(1e-6) == 2 );
    REQUIRE( Approx(formf(2,1,0,1)).epsilon(1e-6) == 1.5 );
    REQUIRE( Approx(formf(2,0,1,1)).epsilon(1e-6) == 1.5 );
    REQUIRE( Approx(formf(2,1,1,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(2,1,2,3)).epsilon(1e-6) == 1.5 );
    REQUIRE( Approx(formf(2,1,0,3)).epsilon(1e-6) == 1.5 );
    REQUIRE( Approx(formf(2,1,4,0)).epsilon(1e-6) == 2 );
    REQUIRE( Approx(formf(2,0,4,8)).epsilon(1e-6) == 0.5 );
  } // GIVEN

  GIVEN( "beryllium oxide input" ){
    REQUIRE( Approx(formf(3,0,0,0)).epsilon(1e-6) == 46.18 );
    REQUIRE( Approx(formf(3,1,0,0)).epsilon(1e-6) == 11.545 );
    REQUIRE( Approx(formf(3,0,1,0)).epsilon(1e-6) == 11.545 );
    REQUIRE( Approx(formf(3,0,0,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(3,1,1,0)).epsilon(1e-6) == 46.18 );
    REQUIRE( Approx(formf(3,1,0,1)).epsilon(1e-6) == 5.67393345 );
    REQUIRE( Approx(formf(3,0,1,1)).epsilon(1e-6) == 5.67393345 );
    REQUIRE( Approx(formf(3,1,1,1)).epsilon(1e-6) == 0 );
    REQUIRE( Approx(formf(3,1,2,3)).epsilon(1e-6) == 29.6660665 );
    REQUIRE( Approx(formf(3,1,0,3)).epsilon(1e-6) == 29.6660665 );
    REQUIRE( Approx(formf(3,1,4,0)).epsilon(1e-6) == 46.18 );
    REQUIRE( Approx(formf(3,0,4,8)).epsilon(1e-6) == 11.545 );
  } // GIVEN

  /* WORK ON TRANSITION BeO next, but for now go eat 
  GIVEN( "beryllium oxide input" ){
    std::vector<std::tuple<int,int,int,int>> inputs 
    { {4,0,0,0}, {4,1,0,0}, {4,0,1,0}, {4,0,0,1}, {4,1,1,0}, {4,1,0,1}, 
      {4,0,1,1}, {4,1,1,1}, {4,1,2,3}, {4,1,0,3}, {4,1,4,0}, {4,0,4,8}, 
      {5,0,0,0}, {5,1,0,0}, {5,0,1,0}, {5,0,0,1}, {5,1,1,0}, {5,1,0,1}, 
      {5,0,1,1}, {5,1,1,1}, {5,1,2,3}, {5,1,0,3}, {5,1,4,0}, {5,0,4,8} };
    for (auto in : inputs){
      REQUIRE( 16 == Approx(formf(get<0>(in),get<1>(in),
                                  get<2>(in),get<3>(in))).epsilon(1e-6) );
    }
  } // GIVEN

  */



  GIVEN( "fcc (like aluminum or lead) lattice input" ){
    std::vector<std::tuple<int,int,int,int>> inputs 
    { {4,0,0,0}, {4,1,0,0}, {4,0,1,0}, {4,0,0,1}, {4,1,1,0}, {4,1,0,1}, 
      {4,0,1,1}, {4,1,1,1}, {4,1,2,3}, {4,1,0,3}, {4,1,4,0}, {4,0,4,8}, 
      {5,0,0,0}, {5,1,0,0}, {5,0,1,0}, {5,0,0,1}, {5,1,1,0}, {5,1,0,1}, 
      {5,0,1,1}, {5,1,1,1}, {5,1,2,3}, {5,1,0,3}, {5,1,4,0}, {5,0,4,8} };
    for (auto in : inputs){
      REQUIRE( 16 == Approx(formf(get<0>(in),get<1>(in),
                                  get<2>(in),get<3>(in))).epsilon(1e-6) );
    }
  } // GIVEN



  GIVEN( "bcc (iron) lattice input" ){
   std::vector<std::tuple<int,int,int,int>> inputs 
     { {6,0,0,0}, {6,1,0,0}, {6,0,1,0}, {6,0,0,1}, {6,1,1,0}, {6,1,0,1},
       {6,0,1,1}, {6,1,1,1}, {6,1,2,3}, {6,1,0,3}, {6,1,4,0}, {6,0,4,8} };
   for (auto in : inputs){
     REQUIRE( 4 == Approx(formf(get<0>(in),get<1>(in),
                                get<2>(in),get<3>(in))).epsilon(1e-6) );
   }
  } // GIVEN


} // TEST CASE
