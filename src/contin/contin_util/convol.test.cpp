#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "convol.h"
#include <iostream>

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}

TEST_CASE( "convol" ){
  GIVEN( "two vectors" ){
    std::vector<double> t1 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0};
    std::vector<double> t2 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double delta = 0.03;
    THEN( "the vectors are correctly convolved and result is returned" ){
      std::vector<double> output = convol( t1, t2, delta );
      std::vector<double> correct = {3.8459762, 2.6993367, 1.0195307, 
        0.53364442, 0.37281623, 0.384, 0.624, 1.008, 1.8, 
        2.16, 0.96, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      equal_vec(output, correct);

           
      t1 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21};
      t2 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21, 0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
      delta = 0.5;
      output = convol( t1, t2, delta );
      correct = {1.1974704E-2, 1.3563056E-2, 1.3531928E-2, 1.379648E-2, 
        1.387114E-2, 1.7875E-2, 2.675E-2, 3.1775E-2, 3.0125E-2, 
        2.52E-2, 1.1025E-2, 0., 0., 0., 0., 0., 0., 0.};

      equal_vec(output, correct);

      t1 = {0.41483349, 0.49122347, 0.28793794, 0.19807373, 0.16013178, 
        0.35027406, 0.54943040};
      t2 = {0.41483349, 0.49122347, 0.28793794, 0.19807373, 0.16013178, 
        0.35027406, 0.5494304, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0};
      delta = 0.3481334;
      output = convol( t1, t2, delta );
      correct = {0.24934471, 0.26630576, 0.25941571, 0.24935982, 0.26861886, 
        0.29515324, 0.28458452, 0.23324588, 0.1398470, 9.5883480E-2, 8.8657E-2,
        0.1004980, 5.2546179E-2, 0., 0., 0., 0., 0., 0., 0., 0.};

      equal_vec(output, correct);

    } // THEN
  } // GIVEN
} // TEST CASE

