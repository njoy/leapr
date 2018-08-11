#include "catch.hpp"
#include "discre/discre_util/prepareParams.h"
#include <range/v3/all.hpp>


TEST_CASE( "prepare parameters helper function" ){
  GIVEN( "inputs" ){
    std::vector<double> ar(50, 0.0), dist(50,0.0), dbw(50,0.0),
      energyNorm(50, 0.0), correctExb(5), energy(2), weights(2), beta(5), 
      exb(5,0.0), betan(5,0.0);
    double tev, tsave, weight, bk, sc;

    energy  = { 0.035, 0.05 };
    weights = { 0.2  , 0.8  };
    beta    = { 0.10, 0.15, 0.30, 0.60, 1.20 };

    tev = 1.723477E-2; bk = 8.617385e-5; sc = 1.0;


    correctExb = {0.904837, 0.860708, 0.740818, 0.548812, 0.301194 };

    std::vector<std::tuple<double,double>> oscEnergiesWeights(energy.size());
    for ( size_t i = 0; i < energy.size(); ++i ){
      oscEnergiesWeights[i] = std::make_tuple(energy[i],weights[i]);
    }

    auto out = prepareParams( oscEnergiesWeights, tev, energyNorm, weight, tsave, ar, dist,
      dbw, bk, exb, betan, beta, sc );
    auto ar_dist_dbw_ranges = std::get<0>(out);
    auto betan_range        = std::get<1>(out);
    auto exb_range          = std::get<2>(out);
    auto ar_range           = std::get<3>(out);

    RANGES_FOR( auto t, ranges::view::zip(betan_range,beta,exb_range,correctExb) ){
      REQUIRE( std::get<1>(t) == Approx(std::get<0>(t)).epsilon(1e-6) );
      REQUIRE( std::get<3>(t) == Approx(std::get<2>(t)).epsilon(1e-6) );
    }


    
    //REQUIRE( 8.213274e-2 == ar_range[0] );
    //REQUIRE( 0.1368162   == ar_range[1] );



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

    REQUIRE( 1.0 == Approx(weight).epsilon(1e-6) );
    REQUIRE( 311.9710021 == Approx(tsave).epsilon(1e-6) );

    REQUIRE( exb.size() == correctExb.size() );
    for( size_t i = 0; i < exb.size(); ++i ){ 
      REQUIRE( correctExb[i] == Approx(exb[i]).epsilon(1e-6) ); 
    }

    REQUIRE( beta.size() == betan.size() );
    for( size_t i = 0; i < beta.size(); ++i ){ 
      REQUIRE( beta[i] == Approx(betan[i]).epsilon(1e-6) ); 
    }  // because sc = 1.0 betan doesn't get scaled

  } // GIVEN
} // TEST CASE

