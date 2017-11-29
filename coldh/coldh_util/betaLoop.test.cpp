#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "betaLoop.h"

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


TEST_CASE( "beta loop helper function" ){
  GIVEN( "inputs" ){
    std::vector<double> alpha { 0.1, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> betan { 0.1, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> rdbex { 1.6, 3.3, 6.6, 20, 5.0, 20, 6.6, 3.3, 1.6, 
      0.0, 0.0 };
    std::vector<double> bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0 };
    std::vector<double> sex { 5, 4, 3, 2, 1, 1, 1, 2, 2, 1, 5 };
    double alpha_val = 0.1;
    double x = 0.85, y = 0.35, swe = 0.32, swo = 0.87;
    double wt = 2.3;
    double tbart = 950;
    int itemp = 0, nbx = 10;
    int a = 0, b = 0;
    int law = 2;
    auto sym_sab = populateSymSab( alpha, betan, true );
    auto sym_sab_2 = populateSymSab( alpha, betan, false );
    betaLoop( betan, rdbex, bex, sex, alpha_val, wt, tbart, x, y, swe, swo, 
      itemp, nbx, a, b, law, sym_sab, sym_sab_2 );

    std::vector<double> correct_sym_sab { 1.70316846, 3.37656089, 5.04886254, 
      6.69968895, 7.86359423, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 
      14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0 };
    std::vector<double> correct_sym_sab_2 { 1.70316846, 1.70258735, 3.37159930,
      3.36200407, 1.68666467, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    
    equal_vec_mega_vec( sym_sab, correct_sym_sab );
    equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );


  } // GIVEN
} // TEST CASE
