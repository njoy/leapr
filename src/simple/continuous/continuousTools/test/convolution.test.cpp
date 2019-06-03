#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simple/continuous/continuousTools/convolution.h"
#include "simple/generalTools/testing.h"


TEST_CASE( "convolution with padding" ){
  std::vector<double> beta0, beta, T1, T2, T3, T4, correctT2, correctT3, correctT4;
  GIVEN( "uniform grid" ){
    WHEN( "Small (length of 3) original beta vec" ){
      beta0 = {-1, 0, 1}; 
      beta  = {-3,-2,-1, 0, 1, 2, 3};

      AND_WHEN( "T1 * T1 -> T2" ){
        T1 = {0, 0, 22.1672, 1, 3, 0, 0};
        correctT2 = {0, 491.38475584, 44.3344, 134.0032, 6, 9,0};
        checkVec( convolutionWithPadding(beta0,beta,T1,T1), correctT2 );
      } // AND WHEN
      AND_WHEN( "T1 * T2 -> T3" ){
        T1 = { 0, 0, 1, 2, 3, 0, 0};
        T2 = { 0, 0, 2, 4, 1, 0, 0};
        correctT3 = { 0, 2, 8,15,14, 3, 0};
        checkVec( convolutionWithPadding(beta0,beta,T1,T2), correctT3 );
      } // AND WHEN
      AND_WHEN( "T1 * Tn -> Tn+1" ){
        beta = {-4,-3,-2,-1, 0, 1, 2, 3, 4 };
        T1   = { 0, 0, 0, 1, 2, 3, 0, 0, 0 };
        correctT2 = { 0, 0, 1, 4,10,12, 9, 0, 0 };
        correctT3 = { 0, 1, 6,21,44,63,54,27, 0 };
        correctT4 = { 1, 8,36,104,214,312,324,216,81 };
        T2 = convolutionWithPadding(beta0,beta,T1,T1); checkVec(T2, correctT2);
        T3 = convolutionWithPadding(beta0,beta,T1,T2); checkVec(T3, correctT3);
        T4 = convolutionWithPadding(beta0,beta,T1,T3); checkVec(T4, correctT4);
      } // AND WHEN
    } // WHEN

    WHEN( "Less small (length of 5) original beta vec" ){
      beta0 = {-2, -1, 0, 1, 2};
      beta = {-4,-3,-2,-1, 0, 1, 2, 3, 4};
      AND_WHEN( "T1 * T1 -> T2" ){
        T1 = { 0, 0, 4, 3, 5, 1, 2, 0, 0};
        correctT2 = {16,24,49,38,47,22,21, 4, 4};
        checkVec( convolutionWithPadding(beta0,beta,T1,T1), correctT2 );
      } // AND WHEN
    } // WHEN
  } // GIVEN
  GIVEN( "beta vecs at uniform vs. nearly uniform" ){
    T1 = {0, 0, 1, 2, 3, 4, 5, 0, 0};
    correctT2 = {1, 4,10,20,35,44,46,40,25};
    WHEN( "T1 * T1 -> T2" ){
      AND_WHEN( "beta grid is uniformly spaced" ){
        beta0 = {-2,-1,0,1,2}; 
        beta  = {-4,-3,-2,-1,0,1,2,3,4};
        THEN( "convolution is correct" ){
          checkVec( convolutionWithPadding(beta0,beta,T1,T1), correctT2 );
        } // THEN
      } // AND WHEN 
      AND_WHEN( "beta grid is almost uniformly spaced" ){
        beta0 = { -2, -0.999999, 0, 0.999999, 2 };
        beta  = { -4,-3,-2,-0.999999, 0, 0.999999, 2, 3, 4 };
        THEN( "convolution is nearly equal to the uniformly spaced result" ){
          checkVec( convolutionWithPadding(beta0,beta,T1,T1), correctT2 );
        } // THEN 
        beta0 = { -2, -1.000001, 0, 1.000001, 2 };
        beta  = { -4,-3,-2,-1.000001, 0, 1.000001, 2, 3, 4 };
        THEN( "convolution is nearly equal to the uniformly spaced result" ){
          checkVec( convolutionWithPadding(beta0,beta,T1,T1), correctT2 );
        } // THEN 

      } // AND WHEN 
    } // WHEN 
  } // GIVEN
  /*
  GIVEN( "vectors" ){
    std::vector<double> beta0 { -4, -3, -2, -1, 0, 1, 2, 3, 4 },
                        beta { -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 },
                        T1 = {  0.0, 0.0, 54.5982, 40.1711, 22.1672, 5.43656, 1, 
                                2, 3, 2, 1, 0.0, 0.0 },
                        T2 = {  0.0, 0.0, 54.5982, 40.1711, 22.1672, 5.43656, 1, 
                                2, 3, 2, 1, 0.0, 0.0 };
    std::vector<double> T3 = convolutionWithPadding(beta0,beta,T1,T2);
    std::vector<double> correctT3 {491.38475584,  44.3344, 134.0032, 6.0, 9.0 };
    for (size_t i = 0; i < T3.size(); ++i){
      //std::cout << T3[i] << std::endl;
      //REQUIRE( correctT3[i] == Approx(T3[i]).epsilon(1e-6) );
    }
  } // GIVEN
  */
} // TEST




TEST_CASE( "convolution" ){
  GIVEN( "vectors" ){

    std::vector<double> beta {0, 1, 3, 6},
                        T1   {1, 2, 5, 0},
                        T2   {3, 1, 3, 2};
    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 1525.7249602 == Approx(T3[0]).epsilon(1e-6) );
  } // GIVEN
  GIVEN( "vectors" ){
    std::vector<double> beta {0, 1, 2, 3},
                        T1   {1,2,3,4},
                        T2   {2,3,4,5};
    std::vector<double> T3 = convolve(beta,T1,T2);
    //for ( auto& x : T3){ std::cout << x << "  ";}
    //REQUIRE( 1525.7249602 == Approx(T3[0]).epsilon(1e-6) );
  } // GIVEN

  GIVEN( "vectors" ){

    std::vector<double> beta {0, 2, 3},
                        T1   {1, 2, 3},
                        T2   {3, 2, 1};
    std::vector<double> T3 = convolve(beta,T1,T2);
    REQUIRE( 154.9252839 == Approx(T3[0]).epsilon(1e-6) );
    REQUIRE( 26.33358414 == Approx(T3[1]).epsilon(1e-6) );
    REQUIRE( 14.0        == Approx(T3[2]).epsilon(1e-6) );
  } // GIVEN
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
    
  } // GIVEN
} // TEST



