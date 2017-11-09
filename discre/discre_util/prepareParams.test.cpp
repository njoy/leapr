#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "prepareParams.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}



TEST_CASE( "prepare parameters helper function" ){
  GIVEN( "inputs" ){
    std::vector<double> ar(50, 0.0), dist(50,0.0), dbw(50,0.0),
                        energyNorm(50, 0.0);
    energyNorm[0] = 2.030778478;
    energyNorm[1] = 2.901112112;
    std::vector<double> energy  { 0.035, 0.05 };
    std::vector<double> weights { 0.2  , 0.8  };
    std::vector<double> beta { 0.10, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> exb     (beta.size(), 0.0);
    std::vector<double> betan   (beta.size(), 0.0);
    double tev = 1.723477E-2;
    double tsave, weight, bk = 8.617385e-5, sc = 1.0;

    prepareParams( energy, weights, tev, energyNorm, weight, tsave, ar, dist,
      dbw, bk, exb, betan, beta, sc );

    std::vector<double> correctExb {0.9512294, 0.9277434, 0.8607079, 
      0.7408182, 0.5488116};
    equal( ar[0], 8.213274e-2 );
    equal( ar[1], 0.1368162   );
    for( auto i = 2; i < ar.size(); ++i ){ equal( ar[i], 0.0 ); }

    equal( dist[0], 4.55739924e-3 );
    equal( dist[1], 2.23263430e-2 );
    for( auto i = 2; i < dist.size(); ++i ){ equal( dist[i], 0.0 ); }

    equal( dbw[0], 0.1282379 );
    equal( dbw[1], 0.3078315 );
    for( auto i = 2; i < dbw.size(); ++i ){ equal( dbw[i], 0.0 ); }

    equal( energyNorm[0], 2.030778 );
    equal( energyNorm[1], 2.901112 );
    for( auto i = 2; i < energyNorm.size(); ++i ){ equal(energyNorm[i], 0.0); }

    equal( weight, 1.0 );
    equal( tsave, 311.9710021 );
    equal_vec( exb, correctExb );
    equal_vec( betan, beta ); // because sc = 1.0 betan doesn't get scaled


  } // GIVEN
} // TEST CASE

