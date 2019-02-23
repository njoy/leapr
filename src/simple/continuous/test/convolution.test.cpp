#include "catch.hpp"
#include "simple/continuous/convolution.h"




TEST_CASE( "convolution" ){
  GIVEN( "vectors" ){

    std::vector<double> beta {0, 1, 3, 6},
                        T1   {1, 2, 5, 0},
                        T2   {3, 1, 3, 2};
    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 1525.7249602 == Approx(T3[0]).epsilon(1e-6) );
  } // GIVEN
  GIVEN( "vectors" ){

    std::vector<double> beta {0, 2, 3},
                        T1   {1, 2, 3},
                        T2   {3, 2, 1};
    std::cout << std::endl;
    std::cout << std::endl;
    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 154.9252839 == Approx(T3[0]).epsilon(1e-6) );
    REQUIRE( 26.33358414 == Approx(T3[1]).epsilon(1e-6) );
    REQUIRE( 14.0        == Approx(T3[2]).epsilon(1e-6) );
    
  } // GIVEN

} // TEST



