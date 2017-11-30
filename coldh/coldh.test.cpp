#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "coldh.h"


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



auto populateSymSab( const std::vector<double>& alpha, const std::vector<double>& beta, bool is_normal ){
  std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
    std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));
  int i = 1;
  for ( auto a = 0; a < alpha.size(); ++a ){
    for ( auto b = 0; b < beta.size(); ++b ){
      if ( is_normal ){ sym_sab[a][b][0] = i; }
      else {sym_sab[a][b][0] = 0; }
      i += 1;
    }
  }
  return sym_sab;
}



TEST_CASE( "coldh" ){
  GIVEN( "inputs" ){
    int itemp = 0;
    double temp = 200.0;
    double tev = 1.723477e-2;
    double sc = 1.0;
    int ncold = 1;
    double tbeta = 2.0;
    double trans_weight = 0.3;
    double scaling = 1.0;
    double dka = 0.2;
    int nbeta = 5, lat = 3;
    std::vector<double> tempf {193093.99765};
    std::vector<double> alpha {0.1, 0.2, 0.4, 0.8, 1.6};
    std::vector<double> beta {0.10, 0.15, 0.30, 0.60, 1.2};
    std::vector<double> ska { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13 };
    auto sym_sab = populateSymSab( alpha, beta, true );
    auto sym_sab_2 = populateSymSab( alpha, beta, false );
    std::vector<double> tempr {200.0};

    coldh( itemp, temp, tev, sc, ncold, trans_weight, tbeta, tempf, tempr, 
      scaling, alpha, beta, dka, ska, nbeta, lat, sym_sab, sym_sab_2 );

    std::vector<double> correctSymSab { 1.7113874, 3.3919859, 5.0714083, 
      6.7348452, 8.4218850, 7.6149802, 8.8664244, 10.115396, 11.306691, 
      12.613772, 13.088162, 14.251586, 15.410194, 16.408195, 17.759561, 
      16.788298, 17.794618, 18.792099, 19.406237, 20.903111, 17.308492, 
      18.066429, 18.809744, 18.777265, 20.590641 };
    std::vector<double> correctSymSab2 { 1.7113874, 2.9207714, 3.7584709, 
      3.7033133, 2.5437978, 7.6149802, 7.6447762, 7.5170871, 6.2697456, 
      3.8402507, 13.088162, 12.312841, 11.499779, 9.2237236, 5.4789194, 
      16.788298, 15.438241, 14.144248, 11.221266, 6.6252107, 17.308492, 
      15.813446, 14.417152, 11.526322, 6.8941688 };
    equal_vec_mega_vec( sym_sab, correctSymSab );
    equal_vec_mega_vec( sym_sab_2, correctSymSab2 );


  } // GIVEN
} // TEST CASEE
