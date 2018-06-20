#include "catch.hpp"
#include "contin/contin_util/checkMomentsEigen.h"
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>



void checkSabEigen( const std::vector<double>& correctSab, 
  const Eigen::Tensor<double,3>& sab ){
	
  REQUIRE( sab.dimension(0)*sab.dimension(1)*sab.dimension(2) == correctSab.size() );

  int l = 0;
  for ( size_t i = 0; i < sab.dimension(0); ++i ){
    for ( size_t j = 0; j < sab.dimension(1); ++j ){
      for ( size_t k = 0; k < sab.dimension(2); ++k ){
        REQUIRE( sab(i,j,k) == Approx(correctSab[l]).epsilon(1e-6) );
	l += 1;
      }
    }
  }

}




TEST_CASE( "check moments eigen" ){

  GIVEN( "inputs" ){

    double sc = 0.99186670867058835;
    std::vector<int> maxt ( 1000, 0.0 );
    maxt[0] = 6;

    std::vector<double> alpha { 1.008e-2, 1.5e-2, 2.52e-2, 3.3e-2, 5.0406e-2 },
      beta { 0.0, 6.375e-3, 1.275e-2, 2.55e-2, 3.825e-2, 5.1e-2, 6.575e-2 };
    
    int itemp = 0, ntempr = 1;
    double f0 = 0.23520650571218535;
    double tbeta = 0.444444;
    double arat = 1, tbar = 1.9344846581861184, explim = -250;

    Eigen::Tensor<double,3> ssm( alpha.size(), beta.size(), ntempr );
    ssm.setZero();

    checkMomentsEigen( sc, alpha, beta, maxt, itemp, f0, tbeta, arat, tbar, ssm );

    std::vector<double> correctSab {0.00000000, 3.04230221, 3.03666664, 
      3.00439217, 2.94493755, 2.85993079, 2.73274581, 0.00000000, 2.49419786, 
      2.49242782, 2.47724974, 2.44682038, 2.40170401, 2.33228656, 0.00000000, 
      1.92380592, 1.92426566, 1.91982023, 1.90827115, 1.88974668, 1.85988639, 
      0.00000000, 1.68058049, 1.68153661, 1.67986667, 1.67343655, 1.66230091, 
      1.64367368, 0.00000000, 1.35861931, 1.35989255, 1.36054304, 1.35866398, 
      1.35426586, 1.34606790 };


    checkSabEigen( correctSab, ssm );


    
  } // GIVEN

  GIVEN( "other inputs" ){

    double sc = 0.99186670867058835;
    std::vector<int> maxt ( 1000, 0.0 );
    maxt[0] = 6; maxt[1] = 1; maxt[2] = 2; maxt[3] = 1; maxt[4] = 2;
    maxt[5] = 1; maxt[6] = 2;

    std::vector<double> alpha { 1.008e-2, 1.5e-2, 2.52e-2, 3.3e-2, 5.0406e-2 },
      beta { 0.0, 6.375e-3, 1.275e-2, 2.55e-2, 3.825e-2, 5.1e-2, 6.575e-2 };
    
    int itemp = 0, ntempr = 1;
    double f0 = 0.23520650571218535;
    double tbeta = 0.444444;
    double arat = 1, tbar = 1.9344846581861184, explim = -250;

    Eigen::Tensor<double,3> ssm( alpha.size(), beta.size(), ntempr );
    ssm.setZero();

    checkMomentsEigen( sc, alpha, beta, maxt, itemp, f0, tbeta, arat, tbar, ssm );

    std::vector<double> correctSab { 0.00000000, 3.04230221, 0.00000000, 
      3.00439217, 0.00000000, 2.85993079, 0.00000000, 0.00000000, 2.49419786, 
      2.49242782, 2.47724974, 2.44682038, 2.40170401, 2.33228656, 0.00000000, 
      1.92380592, 1.92426566, 1.91982023, 1.90827115, 1.88974668, 1.85988639, 
      0.00000000, 1.68058049, 1.68153661, 1.67986667, 1.67343655, 1.66230091, 
      1.64367368, 0.00000000, 1.35861931, 1.35989255, 1.36054304, 1.35866398, 
      1.35426586, 1.34606790 };


    checkSabEigen( correctSab, ssm );

    
  } // GIVEN

} // TEST_CASE
