#include "catch.hpp"
#include "simple/continuous/convolution.h"




TEST_CASE( "convolution" ){
  /*
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
    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 154.9252839 == Approx(T3[0]).epsilon(1e-6) );
    REQUIRE( 26.33358414 == Approx(T3[1]).epsilon(1e-6) );
    REQUIRE( 14.0        == Approx(T3[2]).epsilon(1e-6) );
    
    std::cout << T3[1]*exp(2) << std::endl;
  } // GIVEN
  */
  GIVEN( "vectors" ){

    std::vector<double> beta {0, 2},
                        T1   {1, 3};
    std::vector<double> T2 = convolve(beta,T1,T1);
    REQUIRE( 491.38335029 == Approx(T2[2]*exp(4)).epsilon(1e-6) ); // T2(-4)
    REQUIRE( 66.501504890 == Approx(T2[1]*exp(2)).epsilon(1e-6) ); // T2(-2)
    REQUIRE( 135.00300978 == Approx(T2[0]).epsilon(1e-6) );        // T2(0)
    REQUIRE( 9.0          == Approx(T2[1]).epsilon(1e-6) );        // T2(2)
    REQUIRE( 9.0          == Approx(T2[2]).epsilon(1e-6) );        // T2(4)
    REQUIRE( beta.size()  == 3 );

    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 3931.0668023 == Approx(T3[2]*exp(4)).epsilon(1e-6) ); // T3(-4)
    REQUIRE( 9066.5719883 == Approx(T3[1]*exp(2)).epsilon(1e-6) ); // T3(-2)
    REQUIRE( 21785.154848 == Approx(T3[3]*exp(6)).epsilon(1e-6) ); //T3(-6)
    REQUIRE( 1068.0240782 == Approx(T3[0]).epsilon(1e-6) ); // T3(0)
    REQUIRE( 1227.0270880 == Approx(T3[1]).epsilon(1e-6) ); // T3(2)
    REQUIRE( 72           == Approx(T3[2]).epsilon(1e-6) ); // T3(4)
    REQUIRE( 54           == Approx(T3[3]).epsilon(1e-6) ); // T3(6)
    //std::cout << T3[2] << std::endl;
    //for ( auto x : beta ){ std::cout << x << "  "; }
    //std::cout << std::endl;

    
  } // GIVEN



} // TEST



