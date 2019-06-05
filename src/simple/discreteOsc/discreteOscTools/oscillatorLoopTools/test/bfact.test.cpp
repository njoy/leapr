#define CATCH_CONFIG_MAIN
#include <iostream>
#include "catch.hpp"
#include "simple/discreteOsc/discreteOscTools/oscillatorLoopTools/bfact.h"
#include "generalTools/testing.h"

TEST_CASE( "bfact" ){
  double x, dwc, beta_i, bzero;
  std::vector<double> bPlus (50, 0.0), bMinus (50, 0.0);
  GIVEN( "inputs (small values, reasonable magnitude)" ){
    x = 8.213274e-3; dwc = 1.282379e-2; beta_i = 2.030778478;

    bzero = bfact(x, dwc, beta_i, bPlus, bMinus );

    std::vector<double> nonZeroTerms_minus {1.119176E-2, 6.343540E-5, 2.397033E-7, 
      6.793257E-10, 1.540182E-12, 2.909946E-15, 4.712496E-18, 6.677675E-21, 
      8.410987E-24};
    std::vector<double> nonZeroTerms_plus {1.468732E-3, 1.092496E-6, 5.417594E-10, 
      2.014904E-13, 5.995048E-17, 1.486447E-20, 3.159074E-24, 5.874598E-28};

    THEN( "the outputs bPlus, bMinus, and bzero are correct" ){
      checkVec( bMinus, nonZeroTerms_minus, nonZeroTerms_minus.size() );
      checkVec( bPlus,  nonZeroTerms_plus,  nonZeroTerms_plus.size()  );
      restAreZero( nonZeroTerms_minus.size(), bMinus );
      restAreZero( nonZeroTerms_plus.size(),  bMinus );
      REQUIRE( 0.98727473 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN
  /*
  GIVEN( "inputs that are of medium range size" ){
    double x = 0.01, dwc = 0.2, beta_i = 3.0;
    std::vector<double> bPlus (50, 0.0), bMinus (50, 0.0);
    double bzero = bfact(x, dwc, beta_i, bPlus, bMinus );
    
    std::vector<double>nonZeroTerms_plus = {9.134290E-4, 5.095317E-7, 1.894861E-10, 
      5.285001E-14, 1.179242E-17, 2.192702E-21, 3.494699E-25, 4.873578E-29};
    std::vector<double> nonZeroTerms_minus = {1.834671E-2, 2.055597E-4, 1.535421E-6, 
      8.601593E-9, 3.854963E-11, 1.439728E-13, 4.608864E-16, 1.290968E-18, 
      3.214286E-21, 7.202715E-24};
    THEN( "the outputs bPlus, bMinus, and bzero are correct" ){
      checkVec( bMinus, nonZeroTerms_minus, nonZeroTerms_minus.size() );
      checkVec( bPlus,  nonZeroTerms_plus,  nonZeroTerms_plus.size()  );
      restAreZero( nonZeroTerms_minus.size(), bMinus );
      restAreZero( nonZeroTerms_plus.size(),  bMinus );
      REQUIRE( 0.818751219 == Approx(bzero).epsilon(1e-6) );

    } // THEN
  } // GIVEN


  GIVEN( "some large inputs" ){
    double x = 0.01, dwc = 5.0e+2, beta_i = 3.0e+2;
    double bzero = bfact(x, dwc, beta_i, bPlus, bMinus );
    std::vector<double> bPlus (50, 0.0), bMinus (50, 0.0);
    THEN( "the outputs bPlus, bMinus, and bzero are correct" ){
      restAreZero( bMinus.size(), bMinus );
      restAreZero( bPlus.size(),  bMinus );
      REQUIRE( 0.0 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN

  GIVEN( "large x that invokes the latter case for applying exp terms" ){
  double x = 10.01; dwc = 5.0e+3; beta_i = 3.0;
  double bzero = bfact(x, dwc, beta_i, bPlus, bMinus );
  std::vector<double> bPlus (50, 0.0), bMinus (50, 0.0);
    THEN( "the outputs bPlus, bMinus, and bzero are correct" ){
      restAreZero( bMinus.size(), bMinus );
      restAreZero( bPlus.size(),  bMinus );
      REQUIRE( 0.0 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN
  */
} // TEST CASE
