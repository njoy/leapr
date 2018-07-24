#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "discre/discre_util/oscLoopFuncs_util/bfact.h"
#include "testFuncs.h"
#include <range/v3/all.hpp>


TEST_CASE( "bfact" ){
  double x, dwc, beta_i, bzero;
  std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
  GIVEN( "inputs (small values, reasonable magnitude)" ){
    x = 8.213274e-3; dwc = 1.282379e-2; beta_i = 2.030778478;
    bzero = bfact(x, dwc, beta_i, bplus, bminus );

    double bminusCorrect[9] {1.119176E-2, 6.343540E-5, 2.397033E-7, 
      6.793257E-10, 1.540182E-12, 2.909946E-15, 4.712496E-18, 6.677675E-21, 
      8.410987E-24};
    double bplusCorrect[8] {1.468732E-3, 1.092496E-6, 5.417594E-10, 
      2.014904E-13, 5.995048E-17, 1.486447E-20, 3.159074E-24, 5.874598E-28};

    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      checkAgainstArray( bminusCorrect, bminus );
      checkAgainstArray( bplusCorrect,  bplus  );
      REQUIRE( 0.98727473 == Approx(bzero).epsilon(1e-6) );

    } // THEN
  } // GIVEN
  GIVEN( "inputs that are of medium range size" ){
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    x = 0.01, dwc = 0.2, beta_i = 3.0;
    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    
    double bplusCorrect[8] = {9.134290E-4, 5.095317E-7, 1.894861E-10, 
      5.285001E-14, 1.179242E-17, 2.192702E-21, 3.494699E-25, 4.873578E-29};
    double bminusCorrect[10] = {1.834671E-2, 2.055597E-4, 1.535421E-6, 
      8.601593E-9, 3.854963E-11, 1.439728E-13, 4.608864E-16, 1.290968E-18, 
      3.214286E-21, 7.202715E-24};
    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      checkAgainstArray( bminusCorrect, bminus );
      checkAgainstArray( bplusCorrect,  bplus  );
      REQUIRE( 0.818751219 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN


  auto zero = ranges::view::iota(1,51) 
            | ranges::view::transform([](auto){return 0.0; } );
  GIVEN( "some large inputs" ){
    x = 0.01, dwc = 5.0e+2, beta_i = 3.0e+2;
    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      REQUIRE(ranges::equal(bplus, zero));
      REQUIRE(ranges::equal(bminus,zero));
      REQUIRE( 7.124754E-218 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN

  GIVEN( "large x that invokes the latter case for applying exp terms" ){
    x = 10.01; dwc = 5.0e+3; beta_i = 3.0;
    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      REQUIRE( ranges::equal(bplus, zero) );
      REQUIRE( ranges::equal(bminus,zero) );
      REQUIRE( 0.0 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN
} // TEST CASE
