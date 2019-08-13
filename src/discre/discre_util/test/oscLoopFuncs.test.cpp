#include "catch.hpp"
#include "discre/discre_util/oscLoopFuncs.h"
#include "generalTools/testing.h"
#include <iostream>


TEST_CASE( "negative n terms loop" ){
  int maxdd = 500;
  std::vector<double> bminus (50, 0), bes (maxdd, 0), wts (maxdd,0),
                      wtn (maxdd, 0), ben (maxdd, 0);
 GIVEN( "inputs" ){
    int n = 0, nn = 1;
    double normEnergy = 2.03077847;
    wts[0] = 0.9872747318, wtn[0] = 1.0;

    std::vector<double> bminusStart { 1.119176E-2, 6.343540E-5, 2.397033E-7, 
    6.79325E-10, 1.54018E-12, 2.909946E-15, 4.71249E-18, 6.67767E-21, 
    8.410987E-24 };
    for (size_t i = 0; i < bminusStart.size(); ++i){bminus[i] = bminusStart[i];}

    posNegTerms( n, normEnergy, bminus, wts, wtn, bes, ben, nn, -1 );
    
    THEN( "values are correct" ){
      std::vector<double> 
        besCorrect{0.0, -2.030778, -4.061556, -6.092335},
        wtsCorrect{0.9872747, 1.11917E-2, 6.34354E-5, 2.39703E-7};
      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
      REQUIRE( n == 3 );
    } // THEN
  } // GIVEN

  GIVEN( "inputs" ){
    int n = 4, nn = 8;
    double normEnergy = 0.5; wtn[0] = 0.5;

    std::vector<double> bminusStart { 1e-1, 2e-2, 3e-3, 4e-4, 5e-5, 6e-6, 7e-7, 8e-8, 9e-9 };
    for (size_t i = 0; i < bminusStart.size(); ++i){bminus[i] = bminusStart[i];}

    posNegTerms( n, normEnergy, bminus, wts, wtn, bes, ben, nn, -1 );
    
    THEN( "values are correct" ){
      std::vector<double> 
        besCorrect {0,0,0,0,0,-.5,-1,-1.5,-2,-2.5,-3,-3.5,-4},
        wtsCorrect {0,0,0,0,0,.05,.01,.0015,2E-4,2.5E-5,3E-6,3.5E-7,4E-8};
      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
      REQUIRE( n == 12 );
    } // THEN
  } // GIVEN

  GIVEN( "other inputs" ){
    int n = 0, nn = 1;
    double normEnergy = 2.5;
    wts[0]    = 0.8; wtn[0]    = 1.0;
    bminus[0] = 1.0; bminus[1] = 2.0;
    for ( size_t i = 2; i < 10; ++i ){ bminus[i] = bminus[i-1] + bminus[i-2];}
    
    std::vector<double> 
      besCorrect {0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5, -20, -22.5, -25},
      wtsCorrect {0.8, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 34.0, 55.0, 89.0};
    
    posNegTerms( n, normEnergy, bminus, wts, wtn, bes, ben, nn, -1 );

    THEN( "values are correct" ){
      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
    } // THEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "positive terms" ){
  GIVEN( "inputs" ){
    int maxdd = 500, i = 0, n = 10, nn = 1;
    std::vector<double> bplus (50, 0.0),  bes (maxdd, 0.0), wts (maxdd, 0.0),
                        wtn (maxdd, 0.0), ben (maxdd, 0.0);
    double normalizedEnergy = 2.5;
    wtn[0] = 1.0;
    wts[0] = 0.8; wts[1] = 1.0; wts[2] = 2;

    for ( auto i = 3; i < 11; ++i ){ wts[i] = wts[i-1] + wts[i-2]; }
    for ( auto i = 0; i < 11; ++i ){ bes[i] = -i * 2.5; }

    bplus[0] = 0.08;   bplus[1] = 0.06;  bplus[2] = 0.04;  bplus[3] = 0.02;
    bplus[4] = 0.008;  bplus[5] = 0.006; bplus[6] = 0.004; bplus[7] = 0.002;
    bplus[8] = 0.0008; bplus[9] = 0.0006;
    posNegTerms( n, normalizedEnergy, bplus, wts, wtn, bes, ben, nn, 1 ); 
    
    THEN( "values are correct" ){
      std::vector<double> 
        besCorrect {0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5, -20, -22.5, -25, 
                 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25},
        wtsCorrect {0.8, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 0.08, 0.06, 0.04, 
                 0.02, 8e-3, 6e-3, 4e-3, 2e-3, 8e-4, 6e-4};
      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
     
      REQUIRE( 20 == n );

    } // THEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "oscillator loop" ){
  std::vector<double> wtsCorrect(24), besCorrect(24), wts(24), bes(24);
  double alpha; 

  GIVEN( "inputs" ){
    int maxdd = 500;
    double tbart = 405.894676, temp = 200.0;

    std::vector<double> wts (maxdd, 0.0), bes(maxdd, 0.0), energyNorm(50, 0.0), 
      dbw(50, 0.0), ar(50, 0.0), dist(50,0.0);

    energyNorm[0] = 2.030778;    energyNorm[1] = 2.901112;
    dist[0]       = 4.55739E-3;  dist[1]       = 2.232634E-2;
    dbw[0]        = 0.128237;    dbw[1]        = 0.307831;
    ar[0]         = 8.213274E-2; ar[1]         = 0.1368162;

    alpha  = 0.1;
    std::vector<std::tuple<double,double>> energiesWgts {{0.1,0.2},{0.3,0.8}};

    oscillatorLoop( alpha, dbw, ar, wts, bes, energyNorm, 
      tbart, dist, temp );

    THEN( "ouput is correct" ){
      wtsCorrect = { 0.9573911, 1.085300E-2, 6.151529E-5, 2.324478E-7, 
        1.424275E-3, 1.0594276E-6, 2.793543E-2, 3.166766E-4, 1.794936E-6, 
        4.155853E-5, 3.091273E-8, 4.075663E-4, 4.620178E-6, 2.618736E-8, 
        6.063216E-7, 3.96416E-6, 4.49378E-8, 2.891787E-8, 1.535389E-3, 
        1.740520E-5, 9.865342E-8, 2.284142E-6, 1.231187E-6, 1.395677E-8 };
      besCorrect = { 0.0, -2.030778, -4.061556, -6.092335, 2.030778, 4.061556, 
        -2.901112, -4.931890, -6.962669, -0.8703336, 1.160444, -5.802224, 
        -7.833002, -9.863781, -3.771445, -8.703336, -10.73411, -11.60444, 
        2.901112, 0.8703336, -1.160444, 4.931890, 5.802224, 3.771445 };

      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
  
      REQUIRE( 407.4545311 == Approx(tbart).epsilon(1e-6) );

    } // THEN
  } // GIVEN

  GIVEN( "inputs2" ){
    int maxdd = 500;
    double tbart = 1.8, temp = 296.0;

    std::vector<double> wts (maxdd, 0.0), bes(maxdd, 0.0), energyNorm(50, 0.0), 
      dbw(50, 0.0), ar(50, 0.0), dist(50,0.0);

    energyNorm[0] = 8.0;    energyNorm[1] = 18.0;
    dist[0]       = 0.017;  dist[1]       = 0.079;
    dbw[0]        = 0.021;  dbw[1]        = 0.017;
    ar[0]         = 7e-4;   ar[1]         = 2e-6;

    alpha  = 0.01;

    oscillatorLoop( alpha, dbw, ar, wts, bes, energyNorm, 
      tbart, dist, temp );

    THEN( "ouput is correct" ){
      wtsCorrect = { 0.99962007, 1.9102091E-4, 1.8251428E-8, 6.408037E-8, 
        8.10000E-5, 1.547858E-8 };
      besCorrect = { 0.0, -8.0, -16.0, 8.0, -18.0, -26.0 };

      checkVec(bes,besCorrect,besCorrect.size());
      checkVec(wts,wtsCorrect,wtsCorrect.size());
      restAreZero(besCorrect.size(),bes);
      restAreZero(wtsCorrect.size(),wts);
  
      REQUIRE( 5.5636289 == Approx(tbart).epsilon(1e-6) );

    } // THEN
  } // GIVEN
} // TEST CASE
















