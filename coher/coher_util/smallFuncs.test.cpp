#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "smallFuncs.h"

void equal( double a, double b ){
	std::cout << a << "    " << b << std::endl;
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

void equalHex( std::tuple<int,int,int,double,double>& a, double& b ){
  equal( tausq(std::get<0>(a),std::get<1>(a),std::get<2>(a),std::get<3>(a),std::get<4>(a) ), b );
}

void equalFcc( std::tuple<int,int,int,double>& a, double& b ){
  equal( taufcc(std::get<0>(a),std::get<1>(a),std::get<2>(a),std::get<3>(a) ), b );
}




TEST_CASE( "tausq" ){
  GIVEN( "inputs" ){
    std::vector<std::tuple<int,int,int,double,double>> inputs
      { {0,0,0,2,4}, {1,0,0,2,4}, {0,1,0,2,4}, {0,0,1,2,4}, {1,1,0,2,4}, 
	{1,0,1,2,4}, {0,1,1,2,4}, {1,1,1,2,4}, {1,2,3,4,5}, {5,3,6,4,5},
        {8,7,9,.1,.2} };
    std::vector<double> output { 0, 78.956835, 78.956835, 157.91367041, 
	236.870505, 236.870505, 236.870505, 394.78417604, 2881.924485, 
	14843.885019, 1306.735642 };
    
    for ( auto i = 0; i < output.size(); ++i ){
      equalHex( inputs[i], output[i] );
    }
  } // GIVEN
} // TEST CASE


TEST_CASE( "taufcc" ){
  GIVEN( "inputs" ){
    std::vector<std::tuple<int,int,int,double>> inputs
      { {0,0,0,2}, {1,0,0,2}, {0,1,0,2}, {0,0,1,2}, {1,1,0,2}, 
	{1,0,1,2}, {0,1,1,2}, {1,1,1,2}, {1,2,3,4}, {5,3,6,4},
        {8,7,9,.1} };
    std::vector<double> output { 0, 78.956835, 78.956835, 78.956835,
	210.55156055 };
    
    for ( auto i = 0; i < output.size(); ++i ){
      equalFcc( inputs[i], output[i] );
    }
    
    
  } // GIVEN
} // TEST CASE



    









