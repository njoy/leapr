#include "catch.hpp"
#include "discre/discre_util/prepareParams.h"


TEST_CASE( "prepare parameters helper function" ){
  GIVEN( "inputs" ){
    std::vector<double> ar(50, 0.0), dist(50,0.0), dbw(50,0.0),
      energyNorm(50, 0.0), correctExb(5), energy(2), weights(2), beta(5), 
      exb(5,0.0);
    double tev, tsave, weight, sc;

    energyNorm[0] = 2.030778478;
    energyNorm[1] = 2.901112112;

    tev = 1.723477E-2; sc = 1.0;

    energy  = { 0.035, 0.05 };
    weights = { 0.2  , 0.8  };
    beta    = { 0.10/sc, 0.15/sc, 0.30/sc, 0.60/sc, 1.20/sc };



    std::vector<std::tuple<double,double>> oscEnergiesWeights(energy.size());
    for ( size_t i = 0; i < energy.size(); ++i ){
      oscEnergiesWeights[i] = std::make_tuple(energy[i],weights[i]);
    }
    auto oscEnergiesWeights2 = ranges::view::zip(energy,weights);

    prepareParams( oscEnergiesWeights2, tev, energyNorm, ar, dist,
      dbw, exb, beta );

    correctExb = {0.904837, 0.860708, 0.740818, 0.548812, 0.301194 };

    REQUIRE( 8.213274e-2 == Approx(ar[0]).epsilon(1e-6) );
    REQUIRE( 0.1368162   == Approx(ar[1]).epsilon(1e-6) );

    for( size_t i = 2; i < ar.size(); ++i ){ 
      REQUIRE( 0.0 == Approx(ar[i]).epsilon(1e-6) ); 
    }

    REQUIRE( 4.55739924e-3 == Approx(dist[0]).epsilon(1e-6) );
    REQUIRE( 2.23263430e-2 == Approx(dist[1]).epsilon(1e-6) );
    for( size_t i = 2; i < dist.size(); ++i ){ 
      REQUIRE( 0.0 == Approx(dist[i]).epsilon(1e-6) ); 
    }


    REQUIRE( 0.1282379 == Approx(dbw[0]).epsilon(1e-6) );
    REQUIRE( 0.3078315 == Approx(dbw[1]).epsilon(1e-6) );
    for( size_t i = 2; i < dbw.size(); ++i ){ 
      REQUIRE( 0.0 == Approx(dbw[i]).epsilon(1e-6) ); 
    }


    REQUIRE( 2.030778 == Approx(energyNorm[0]).epsilon(1e-6) );
    REQUIRE( 2.901112 == Approx(energyNorm[1]).epsilon(1e-6) );
    for( size_t i = 2; i < energyNorm.size(); ++i ){ 
      REQUIRE( 0.0 == Approx(energyNorm[i]).epsilon(1e-6) ); 
    }

    //REQUIRE( 1.0 == Approx(weight).epsilon(1e-6) );
    //REQUIRE( 311.9710021 == Approx(tsave).epsilon(1e-6) );

    REQUIRE( exb.size() == correctExb.size() );
    for( size_t i = 0; i < exb.size(); ++i ){ 
      REQUIRE( correctExb[i] == Approx(exb[i]).epsilon(1e-6) ); 
    }

  } // GIVEN
} // TEST CASE

