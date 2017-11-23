#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "coldh.h"

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
    auto sym_sab = populateSymSab( alpha, beta );
    std::vector<double> tempr {200.0};

    coldh( itemp, temp, tev, sc, ncold, trans_weight, tbeta, tempf, tempr, 
      scaling, alpha, beta, dka, ska, nbeta, lat, sym_sab );


  } // GIVEN
} // TEST CASE
