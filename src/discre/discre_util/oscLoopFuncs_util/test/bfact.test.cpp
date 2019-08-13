#define CATCH_CONFIG_MAIN
#include <iostream>
#include "catch.hpp"
#include "discre/discre_util/oscLoopFuncs_util/bfact.h"
#include "generalTools/testing.h"
#include <range/v3/all.hpp>


TEST_CASE( "bfact" ){
  double x, dwc, beta_i, bzero;
  std::vector<double> bplus (50, 0.0), bminus (50, 0.0);

  GIVEN( "inputs (small values, reasonable magnitude)" ){
    x = 8.213274e-3; dwc = 1.282379e-2; beta_i = 2.030778478;

    bzero = bfact(x, dwc, beta_i, bplus, bminus );

    std::vector<double> 
      bminusStart {1.11917E-2, 6.34354E-5, 2.39703E-7, 6.79325E-10, 1.54018E-12, 
      2.909946E-15, 4.712496E-18, 6.677675E-21, 8.410987E-24},
      bplusStart {1.46873E-3, 1.09249E-6, 5.41759E-10, 2.01490E-13, 5.99504E-17, 
      1.48644E-20, 3.15907E-24, 5.87459E-28};

    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      checkVec(bminus,bminusStart,bminusStart.size());
      restAreZero(bminusStart.size(),bminus);
      checkVec(bplus,bplusStart,bplusStart.size());
      restAreZero(bplusStart.size(),bplus);
      REQUIRE( 0.98727473 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN

  GIVEN( "inputs that are of medium range size" ){
    double x = 0.01, dwc = 0.2, beta_i = 3.0;
    std::vector<double> bplus(50,0.0), bminus(50,0.0);
    double bzero = bfact(x, dwc, beta_i, bplus, bminus );
    
    std::vector<double> 
      bminusStart = {1.83467E-2, 2.05559E-4, 1.53542E-6, 8.60159E-9, 3.85496E-11, 
      1.439728E-13, 4.60886E-16, 1.29096E-18, 3.21428E-21, 7.20271E-24},
      bplusStart = {9.134290E-4, 5.095317E-7, 1.894861E-10, 5.285001E-14, 
      1.179242E-17, 2.192702E-21, 3.494699E-25, 4.873578E-29};

    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      checkVec(bminus,bminusStart,bminusStart.size());
      restAreZero(bminusStart.size(),bminus);
      checkVec(bplus,bplusStart,bplusStart.size());
      restAreZero(bplusStart.size(),bplus);
      REQUIRE( 0.818751219 == Approx(bzero).epsilon(1e-6) );
    } // THEN
  } // GIVEN

  GIVEN( "large x that invokes the latter case for applying exp terms" ){
    double x = 10.01; dwc = 5.0e+3; beta_i = 3.0;
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    double bzero = bfact(x, dwc, beta_i, bplus, bminus );
    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      restAreZero(0,bminus);
      restAreZero(0,bplus);
      REQUIRE( 0.0 == Approx(bzero).epsilon(1e-6) );
    } // THEN

    x = 3.76; dwc = 0.03; beta_i = 2.0;
    for (auto& x : bplus ){ x = 0; } 
    for (auto& x : bminus){ x = 0; }
    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    std::vector<double> 
      bminusStart = { 20.709925, 36.005053, 48.907781, 53.896867, 49.665579, 
      39.19101, 26.98526, 16.45973, 9.003310, 4.461336, 2.019739, 0.841428, 
      0.324569, 0.116543, 3.91362E-2, 1.23414E-2, 3.66803E-3, 1.03083E-3, 
      2.74726E-4, 6.96161E-5, 1.68134E-5, 3.87866E-6, 8.56348E-7, 1.81281E-7, 
      3.68572E-8, 7.20827E-9, 1.35801E-9, 2.46784E-10, 4.33130E-11, 7.35035E-12, 
      1.20742E-12, 1.92183E-13, 2.96682E-14 },
      bplusStart = { 2.802783690, 0.659455550, 0.121230269, 1.808038E-2, 
      2.25481E-3, 2.40797E-4, 2.24390E-5, 1.85229E-6, 1.37120E-7, 9.19550E-9, 
      5.63399E-10, 3.17650E-11, 1.65825E-12, 8.05826E-14, 3.66221E-15, 
      1.56294E-16, 6.28667E-18, 2.39104E-19, 8.62401E-21, 2.95754E-22, 
      9.66691E-24, 3.01803E-25, 9.01787E-27, 2.58356E-28, 7.10884E-30 };

    THEN( "the outputs bplus, bminus, and bzero are correct" ){
      checkVec(bminus,bminusStart,bminusStart.size());
      restAreZero(bminusStart.size(),bminus);
      checkVec(bplus,bplusStart,bplusStart.size());
      restAreZero(bplusStart.size(),bplus);

      REQUIRE( 8.925283  == Approx(bzero).epsilon(1e-6) );
    } // THEN


  } // GIVEN
} // TEST CASE
