#define CATCH_CONFIG_MAIN
#include "../catch.hpp" 
#include "skold.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-5 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-5 );
}

void equal_vec_mega_vec( std::vector<std::vector<std::vector<double>>> a, 
  std::vector<double> b ){
  REQUIRE( a.size()*a[0].size()*a[0][0].size() == b.size() );
  int i = 0;
  for ( auto a1 : a ){
    for ( auto a2 : a1 ){
      for ( auto a3 : a2 ){
        equal( a3, b[i] );
        i += 1;
      }
    }
  }
}

auto populateSymSab( int a1, int b1 ){
  std::vector<std::vector<std::vector<double>>> sym_sab(a1,
    std::vector<std::vector<double>>(b1,std::vector<double>(1,0.0)));
  int i = 1;
  for ( auto a = 0; a < a1; ++a ){
    for ( auto b = 0; b < b1; ++b ){
      sym_sab[a][b][0] = i;
      i += 1;
    }
  }
  return sym_sab;
}



TEST_CASE( "skold" ){
    double cfrac = 0.3;
    auto symsab = populateSymSab( 5, 5 );
    double temp = 200.0, awr = 8.9;
    int nalpha = 5, nbeta = 5, ntempr = 1, itemp = 0;
    int lat = 3, nka = 3;
    double dka = 2.3, scaling = 1.0;
    std::vector<double> alpha  { 0.10, 0.20, 0.40, 0.80, 1.60 };
    std::vector<double> beta   { 0.10, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> skappa { 1.10, 2.20, 3.30, 4.40, 5.50 };
    std::vector<double> rightSymSab { 0.7748423, 1.695805, 2.725674, 
      3.834036, 5, 4.525767, 5.779165, 7.133967, 8.547886, 10, 9.056398, 
      10.95454, 12.81736, 14.66346, 16.5, 16.62805, 19.17417, 21.51546, 
      23.77780, 26, 21, 22, 23, 24, 25 };
  GIVEN( "inputs" ){
    skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, dka, scaling, symsab );
    equal_vec_mega_vec( symsab, rightSymSab );
  } // GIVEN

  GIVEN( "small values in skappa" ){
    skappa = { 0.001, 0.0012, 0.0013 };
    WHEN( "value of dka is small" ){
      THEN( "no change to scattering law" ){
        for ( int i = 0; i < rightSymSab.size(); ++i ){ rightSymSab[i] = i+1; }
        dka = 0.03;
        skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka,
          dka, scaling, symsab );
        equal_vec_mega_vec( symsab, rightSymSab );
        dka = 1.03;
        skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, 
          dka, scaling, symsab );
        equal_vec_mega_vec( symsab, rightSymSab );
      } // THEN
    } // WHEN

    WHEN( "value of dka is of moderate size" ){
      THEN( "the scattering law is correctly amended" ){

        dka = 2.03;
        rightSymSab = { 0.7362570, 1.435067, 2.134116, 2.833354, 3.532743, 
          4.248877, 4.946653, 5.644851, 6.343379, 7.042169, 11, 12, 13, 14, 15, 
          16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
        skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, 
          dka, scaling, symsab );
        equal_vec_mega_vec( symsab, rightSymSab );
        
        symsab = populateSymSab( 5, 5 );
        cfrac = 1.3;
        rightSymSab = { -0.1428862, -0.4480389, -0.7521607, -1.055464, 
          -1.358109, -1.588199, -1.897833, -2.205644, -2.512024, -2.817265, 11,
          12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
        skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, 
          dka, scaling, symsab );
        equal_vec_mega_vec( symsab, rightSymSab );

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
        skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, 
          dka, scaling, symsab );
        equal_vec_mega_vec( symsab, rightSymSab );

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
