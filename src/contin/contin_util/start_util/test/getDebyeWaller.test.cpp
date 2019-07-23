#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin_util/start_util/getDebyeWaller.h"
#include "range/v3/all.hpp"


template <typename Float>
auto makeGrid( int len, Float delta ){
  return ranges::view::iota(0,len) 
       | ranges::view::transform([delta](auto x){ return x*delta;});
}

TEST_CASE( "getting Debye Waller coefficient" ){
  GIVEN( "a specified beta grid" ){ 
    std::vector<double> p1  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                        p2  {0.01, 0.02, 0.03, 0.04, 0.05, 0.06},
                        betaGrid1 = makeGrid(int(p1.size()),2.0),
                        betaGrid2 = makeGrid(int(p2.size()),0.1);
    WHEN( "n = 0 (used for debye-waller coefficient) " ){
      THEN( "returned value has as most a 1e-6 percent error" ){
        auto beta_P_1 = ranges::view::zip(betaGrid1,p1);
        auto beta_P_2 = ranges::view::zip(betaGrid2,p2);
        REQUIRE( getDebyeWaller(beta_P_1) == Approx(1444532.840).epsilon(1e-6) );
        REQUIRE( getDebyeWaller(beta_P_2) == Approx(0.035514341).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
