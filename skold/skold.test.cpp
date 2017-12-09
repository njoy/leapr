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
  GIVEN( "inputs" ){
    double cfrac = 0.3;
    auto symsab = populateSymSab( 5, 5 );
    double temp = 200.0, awr = 8.9;
    int nalpha = 5, nbeta = 5, ntempr = 1, itemp = 0;
    int lat = 3, nka = 3;
    double dka = 0.02, scaling = 1.0;
    std::vector<double> alpha  { 0.10, 0.20, 0.40, 0.80, 1.6 };
    std::vector<double> beta   { 0.10, 0.15, 0.30, 0.60, 1.2 };
    std::vector<double> skappa { 1.20, 2.30, 3.40 };
    skold( cfrac, itemp, temp, alpha, beta, skappa, ntempr, awr, lat, nka, dka, scaling, symsab );
    REQUIRE( true );
  



  } // GIVEN
} // TEST CASE
