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
    std::vector<double> rdbex { 1.6, 3.3, 6.6, 20, 5.0, 20, 6.6, 3.3, 1.6, 0.0, 0.0 };
    std::vector<double> bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 };
    std::vector<double> sex { 5, 4, 3, 2, 1, 1, 1, 2, 2, 1, 5 };
    double alpha_val = 0.1;
    double x = 0.85292696, y = 0.35;
    double wt = 2.3;
    double tbart = 950;
    int itemp = 0, nbx = 10;
    int a = 0, b = 0;
    int law = 2;
    auto sym_sab = populateSymSab( alpha, betan );
    betaLoop( betan, rdbex, bex, sex, alpha_val, wt, tbart, x, y, itemp, nbx, a, b, law, sym_sab );



  } // GIVEN
} // TEST CASE
