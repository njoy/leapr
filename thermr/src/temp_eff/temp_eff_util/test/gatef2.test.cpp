#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "temp_eff/temp_eff_util/gatef2.h"

TEST_CASE( "gatef2" ){
  GIVEN( "two vectors" ){
    std::vector<double> 
      t1 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0},
      t2 {0.2, 0.6, 0.8, 2.0, 6.0, 8.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      output(18), correct(18); 
    double delta = 0.03;

    THEN( "the vectors are correctly convolved and result is returned" ){

      output = convol( t1, t2, delta ),
      correct = {3.8459762, 2.6993367, 1.0195307, 0.53364442, 0.37281623, 
        0.384, 0.624, 1.008, 1.8, 2.16, 0.96, 0, 0, 0, 0, 0, 0, 0};

      REQUIRE( output.size() == correct.size() );
      for ( size_t i = 0; i < output.size(); ++i ){
        REQUIRE( output[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
           
      t1 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21};
      t2 = {0.01, 0.04, 0.09, 0.11, 0.16, 0.21, 0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};

      delta = 0.5;
      output = convol( t1, t2, delta );
      correct = {1.1974704E-2, 1.3563056E-2, 1.3531928E-2, 1.379648E-2, 
        1.387114E-2, 1.7875E-2, 2.675E-2, 3.1775E-2, 3.0125E-2, 
        2.52E-2, 1.1025E-2, 0., 0., 0., 0., 0., 0., 0.};

      REQUIRE( output.size() == correct.size() );
      for ( size_t i = 0; i < output.size(); ++i ){
        REQUIRE( output[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
 
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

      REQUIRE( output.size() == correct.size() );
      for ( size_t i = 0; i < output.size(); ++i ){
        REQUIRE( output[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
      }
 
    } // THEN
  } // GIVEN
} // TEST CASE

