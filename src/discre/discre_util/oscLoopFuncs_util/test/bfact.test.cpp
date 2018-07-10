#define CATCH_CONFIG_MAIN
#include <iostream>
#include "catch.hpp"
#include "discre/discre_util/oscLoopFuncs_util/bfact.h"
#include <range/v3/all.hpp>


template <typename arrayT, typename rangeT1, typename rangeT2> 
auto checkBminusBplus( arrayT bmRef, arrayT bpRef, rangeT1 bm, rangeT2 bp ){

  RANGES_FOR( auto t, ranges::view::zip( bm, bmRef ) ){
    REQUIRE( std::get<1>(t) == Approx(std::get<0>(t)).epsilon(1e-6) );
  }
  RANGES_FOR( auto x, bm | ranges::view::slice(int(bmRef.size()),ranges::end )){ 
    REQUIRE( 0.0 == Approx(x).epsilon(1e-6) ); 
  }
  RANGES_FOR( auto t, ranges::view::zip( bp, bpRef ) ){
    REQUIRE( std::get<1>(t) == Approx(std::get<0>(t)).epsilon(1e-6) );
  }
  RANGES_FOR( auto x, bp | ranges::view::slice(int(bpRef.size()),ranges::end )){ 
    REQUIRE( 0.0 == Approx(x).epsilon(1e-6) ); 
  }

}


TEST_CASE( "bfact" ){
  GIVEN( "vectors as inputs" ){
    double x, dwc, beta_i, bzero;
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    WHEN( "inputs (small values, reasonable magnitude)" ){
      x = 8.213274e-3; dwc = 1.282379e-2; beta_i = 2.030778478;

      //bzero = bfact(x, dwc, beta_i, bplus, bminus );
      auto out = bfact(x, dwc, beta_i, bplus, bminus );
      double bzero  = std::get<0>(out);
      auto b_minus  = std::get<1>(out);
      auto b_plus   = std::get<2>(out);

      //std::cout << b_plus << std::endl;
      //std::cout << b_minus << std::endl;
      std::array<double,9> bminusTrue {{1.119176E-2, 6.343540E-5, 2.397033E-7, 
        6.793257E-10, 1.540182E-12, 2.909946E-15, 4.712496E-18, 6.677675E-21, 
        8.410987E-24}},
      bplusTrue {{1.468732E-3, 1.092496E-6, 5.417594E-10, 2.014904E-13, 
        5.995048E-17, 1.486447E-20, 3.159074E-24, 5.874598E-28, 0.0}};

      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        checkBminusBplus( bminusTrue, bplusTrue, b_minus, b_plus );
        REQUIRE( 0.98727473 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "inputs that are of medium range size" ){
      double x = 0.01, dwc = 0.2, beta_i = 3.0;
      std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      auto out = bfact(x, dwc, beta_i, bplus, bminus );
      double bzero  = std::get<0>(out);
      auto b_minus  = std::get<1>(out);
      auto b_plus   = std::get<2>(out);
    
      std::array<double,10> bplusTrue {{9.134290E-4, 5.095317E-7, 1.894861E-10, 
        5.285E-14, 1.179242E-17, 2.1927E-21, 3.494699E-25, 4.873578E-29, 0, 0}},
      bminusTrue {{1.834671E-2, 2.055597E-4, 1.535421E-6, 8.601593E-9, 
        3.854963E-11, 1.439728E-13, 4.608864E-16, 1.290968E-18, 3.214286E-21, 
        7.202715E-24}};
      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        checkBminusBplus( bminusTrue, bplusTrue, b_minus, b_plus );
        REQUIRE( 0.818751219 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "some large inputs" ){
      double x = 0.01, dwc = 4.0e+2, beta_i = 1.5e+2;
      std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      auto out = bfact(x, dwc, beta_i, bplus, bminus );
      double bzero  = std::get<0>(out);
      auto b_minus  = std::get<1>(out);
      auto b_plus   = std::get<2>(out);


      std::array<double,10> bminusTrue {{0.0, 0.0, 0.0, 0.0, 3.61666671336e-25, 
        112515.6977012, 3.000343793e34, 7.0006283212e63, 1.4519462116e93, 
        2.7102326105e122}},
        bplusTrue {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};

      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        checkBminusBplus( bminusTrue, bplusTrue, b_minus, b_plus );
        REQUIRE( 7.124754E-218 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "large x that invokes the latter case for applying exp terms" ){
      double x = 10.01; dwc = 4.0e+3; beta_i = 3.0;
      std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      auto out = bfact(x, dwc, beta_i, bplus, bminus );
      double bzero  = std::get<0>(out);
      auto b_minus  = std::get<1>(out);
      auto b_plus   = std::get<2>(out);

      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        for(auto entry : b_plus ){REQUIRE(0.0 == Approx(entry).epsilon(1e-6) );}
        for(auto entry : b_minus){REQUIRE(0.0 == Approx(entry).epsilon(1e-6) );}
        REQUIRE( 5.443974920E-171 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN



  GIVEN( "ranges as inputs" ){
    double x, dwc, beta_i, bzero;
    auto bplus = ranges::view::iota(0,50) 
               | ranges::view::transform([](auto){return 0.0; } );
    auto bminus = bplus;

    //std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    WHEN( "inputs (small values, reasonable magnitude)" ){
      x = 8.213274e-3; dwc = 1.282379e-2; beta_i = 2.030778478;

    /*
      bzero = bfact(x, dwc, beta_i, bplus, bminus );

      std::vector<double> bminusTrue = {1.119176E-2, 6.343540E-5, 2.397033E-7, 
        6.793257E-10, 1.540182E-12, 2.909946E-15, 4.712496E-18, 6.677675E-21, 
        8.410987E-24};
      std::vector<double> bplusTrue = {1.468732E-3, 1.092496E-6, 5.417594E-10, 
        2.014904E-13, 5.995048E-17, 1.486447E-20, 3.159074E-24, 5.874598E-28};

      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        checkBminusBplus( bminusTrue, bplusTrue, bminus, bplus );
        REQUIRE( 0.98727473 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "inputs that are of medium range size" ){
      double x = 0.01, dwc = 0.2, beta_i = 3.0;
      std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      double bzero = bfact(x, dwc, beta_i, bplus, bminus );
    
      std::vector<double>bplusTrue = {9.134290E-4, 5.095317E-7, 1.894861E-10, 
        5.285001E-14, 1.179242E-17, 2.192702E-21, 3.494699E-25, 4.873578E-29};
      std::vector<double> bminusTrue = {1.834671E-2, 2.055597E-4, 1.535421E-6, 
        8.601593E-9, 3.854963E-11, 1.439728E-13, 4.608864E-16, 1.290968E-18, 
        3.214286E-21, 7.202715E-24};
      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        checkBminusBplus( bminusTrue, bplusTrue, bminus, bplus );
        REQUIRE( 0.818751219 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN


    WHEN( "some large inputs" ){
      double x = 0.01, dwc = 5.0e+2, beta_i = 3.0e+2;
      double bzero = bfact(x, dwc, beta_i, bplus, bminus );
      std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        for(auto entry : bplus ){REQUIRE( 0.0 == Approx(entry).epsilon(1e-6) );}
        for(auto entry : bminus){REQUIRE( 0.0 == Approx(entry).epsilon(1e-6) );}
        REQUIRE( 7.124754E-218 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "large x that invokes the latter case for applying exp terms" ){
    double x = 10.01; dwc = 5.0e+3; beta_i = 3.0;
    double bzero = bfact(x, dwc, beta_i, bplus, bminus );
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
      THEN( "the outputs bplus, bminus, and bzero are correct" ){
        for(auto entry : bplus ){REQUIRE(0.0 == Approx(entry).epsilon(1e-6) );}
        for(auto entry : bminus){REQUIRE(0.0 == Approx(entry).epsilon(1e-6) );}
        REQUIRE( 0.0 == Approx(bzero).epsilon(1e-6) );
      } // THEN
    */
    } // WHEN
  } // GIVEN

} // TEST CASE
