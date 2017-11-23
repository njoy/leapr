#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "betaLoop.h"

auto populateSymSab( const std::vector<double>& alpha, const std::vector<double>& beta ){
  std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
    std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));
  int i = 1;
  for ( auto a = 0; a < alpha.size(); ++a ){
    for ( auto b = 0; b < beta.size(); ++b ){
      sym_sab[a][b][0] = i;
      i += 1;
    }
  }
  return sym_sab;
}



TEST_CASE( "beta loop helper function" ){
  GIVEN( "inputs" ){
    std::vector<double> alpha { 0.1, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> betan { 0.1, 0.15, 0.30, 0.60, 1.20 };
    double alp = 0.23;
    double x = 0.85292696;
    int itemp = 0;
    int a = 0;
    int law = 2;
    auto sym_sab = populateSymSab( alpha, betan );
    betaLoop( betan, alp, x, itemp, a, law, sym_sab );



  } // GIVEN
} // TEST CASE
