#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "skold/skold.h"
#include "generalTools/testing.h"

TEST_CASE( "skold" ){
  double cfrac = 0.3, tev = 0.0172348, awr = 8.9, dka = 2.3, scaling = 1.0;
  int nalpha = 5, nbeta = 5, itev = 0;
  std::vector<double> alpha { 0.10, 0.20, 0.40, 0.80, 1.60 }, 
    beta { 0.1, 0.15, 0.3, 0.6, 1.2 }, skappa { 1.1, 2.2, 3.3, 4.4, 5.5 }, 
    rightSymSab { 0.7748423, 1.695805, 2.725674, 3.834036, 5, 4.525767, 
    5.779165, 7.133967, 8.547886, 10, 9.056398, 10.95454, 12.81736, 14.66346,
    16.5, 16.62805, 19.17417, 21.51546, 23.77780, 26, 21, 22, 23, 24, 25 };
  std::vector<double> symSab = ranges::view::iota(1,26);

  GIVEN( "medium sized skappa values" ){
    WHEN( "dka is a bit large" ){
      THEN( "the scattering law is correctly changed" ){
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
      } // THEN
    } // WHEN
    WHEN( "dka is really big" ){
      THEN( "the scattering law is correctly changed" ){
        dka = 25.0; cfrac = 1.3;
        rightSymSab = { 0.6489166, 1.614790, 2.691405, 3.827429, 5, 3.545016, 
          5.406339, 7.019986, 8.535378, 10, 11.53434, 13.04755, 14.54618, 
          16.03430, 17.51459, 19.78249, 21.38301, 22.97954, 24.57282, 26.16343,
          29.73225, 31.50051, 33.26723, 35.03263, 36.79689 };
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "small values in skappa" ){
    skappa = { 0.001, 0.0012, 0.0013 };
    WHEN( "value of dka is small" ){
      THEN( "no change to scattering law" ){
        for ( size_t i = 0; i < rightSymSab.size(); ++i ){ rightSymSab[i] = i+1; }
        dka = 0.03;
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
        dka = 1.03;
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
      } // THEN
    } // WHEN

    WHEN( "value of dka is of moderate size" ){
      THEN( "the scattering law is correctly amended" ){

        dka = 2.03;
        rightSymSab = { 0.7362570, 1.435067, 2.134116, 2.833354, 3.532743, 
          4.248877, 4.946653, 5.644851, 6.343379, 7.042169, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );

        symSab = ranges::view::iota(1,26);
        cfrac = 1.3;
        rightSymSab = { -0.1428862, -0.4480389, -0.7521607, -1.055464, 
          -1.358109, -1.588199, -1.897833, -2.205644, -2.512024, -2.817265, 11,
          12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
      } // THEN
    } // WHEN
    WHEN( "value of dka is large" ){
      THEN( "the scattering law is correctly amended" ){
        dka = 5.03;
        cfrac = 1.3;
        rightSymSab = { -0.1528274, -0.4579653, -0.7620900, -1.065412, 
          -1.368090, -1.602126, -1.911557, -2.219208, -2.525467, -2.830620, 
          -3.032806, -3.348865, -3.661900, -3.972597, -4.281455, -4.441460, 
          -4.767483, -5.088515, -5.405735, -5.719990, 21, 22, 23, 24, 25 };
        skold( cfrac, tev, alpha, beta, skappa, awr, dka, scaling, symSab );
        REQUIRE( ranges::equal(symSab, rightSymSab, equal) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
