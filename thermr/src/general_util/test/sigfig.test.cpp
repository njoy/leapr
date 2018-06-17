#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "general_util/sigfig.h"


TEST_CASE( "sigfig" ){
  GIVEN( "ins" ){
    double out;
    std::vector<std::tuple<int,int,double>> in { {1,-3,-2.0}, {1,0,1.0}, 
      {1,1,2.0}, {2,1,1.3}, {3,-1,1.22}, {3,1,1.24}, {3,2,1.25}, {3,-2,1.21}, 
      {5,0,1.2346}, {8,-1,1.2345678} };
    for ( size_t i = 0; i < in.size(); ++i ){
      out = sigfig( 1.234567890123456789, std::get<0>(in[i]), std::get<1>(in[i]) );
      REQUIRE( std::get<2>(in[i]) == Approx(out).epsilon(1e-8) );
    }
  } // GIVEN
  GIVEN( "ins" ){
    double out;
    std::vector<std::tuple<int,int,double>> in { {5,0,2.1212}, {7,0,2.121212}, 
      {7,-5,2.121207}, {10,0,2.12121201} };
    for ( size_t i = 0; i < in.size(); ++i ){
      out = sigfig( 2.121212, std::get<0>(in[i]), std::get<1>(in[i]) );
      REQUIRE( std::get<2>(in[i]) == Approx(out).epsilon(1e-8) );
    }

  } // GIVEN

} // TEST CASE

