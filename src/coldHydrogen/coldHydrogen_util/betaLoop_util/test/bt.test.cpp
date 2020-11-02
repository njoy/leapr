#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coldHydrogen/coldHydrogen_util/betaLoop_util/bt.h"


TEST_CASE( "bt" ){
  int j;
  double pj, x;
  using std::get;

  GIVEN( "odd input value for j" ){
    std::vector<std::tuple<int,double>> inputs { {1,0.85}, {1,0.35}, {5,0.35} };
    std::vector<double> output { 0.48388278, 0.34887661, 9.52577E-3 };
    for (size_t i = 0; i < inputs.size(); ++i){
      auto in = inputs[i];
      REQUIRE( bt(get<0>(in),get<1>(in)) == Approx(output[i]).epsilon(1e-6) );
    }
  } // GIVEN

  GIVEN( "even input value for j" ){
    std::vector<std::tuple<int,double>> inputs { {2,3.85}, {4,3.85}, {6,0.005} };
    std::vector<double> output { 2.40889540E-5, 8.5675066E-17, 4.76258212E-2 };
    for (size_t i = 0; i < inputs.size(); ++i){
      auto in = inputs[i];
      REQUIRE( bt(get<0>(in),get<1>(in)) == Approx(output[i]).epsilon(1e-6) );
    }
  } // GIVEN
} // TEST CASE
