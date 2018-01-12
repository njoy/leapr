#include "catch.hpp"
#include "discre/discre_util/oscLoopFuncs.h"

void equal3( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal3_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( size_t i = 0; i < a.size(); ++i ){
    equal3( a[i], b[i] );
  }
}

TEST_CASE( "negative n terms loop" ){
  int maxdd = 500;
  std::vector<double> bminus (50, 0.0), bes (maxdd, 0.0), wts (maxdd, 0.0),
                      wtn (maxdd, 0.0), ben (maxdd, 0.0);
 GIVEN( "inputs" ){
    int n = 0, nn = 1;
    double wtsn = 0.987274731, besn = 0.0, normEnergy = 2.03077847;
    wts[0] = 0.9872747318, wtn[0] = 1.0;

    bminus[0] = 1.119176E-2; bminus[1] = 6.343540E-5; bminus[2] = 2.397033E-7; 
    bminus[3] = 6.79325E-10; bminus[4] = 1.54018E-12; bminus[5] = 2.909946E-15; 
    bminus[6] = 4.71249E-18; bminus[7] = 6.67767E-21; bminus[8] = 8.410987E-24;

    posNegTerms( n, normEnergy, bminus, wts, wtn, bes, ben, nn, -1 );
    
    THEN( "values are correct" ){
      equal3( n, 3 );

      std::vector<double> besVals {0.0, -2.030778, -4.061556, -6.092335};
      for ( size_t i = 0; i < bes.size(); ++i ){ 
        if ( i < 4 ){ equal3( bes[i], besVals[i] ); }
        else { equal3( bes[i], 0.0 ); }
      }

      std::vector<double> wtsVals {0.9872747, 1.119176E-2, 6.343540E-5, 
        2.397033E-7};
      for ( size_t i = 0; i < wts.size(); ++i ){
        if ( i < 4 ){ equal3( wts[i], wtsVals[i] ); }
        else { equal3( wts[i], 0.0 ); }
      }
    } // THEN
  } // GIVEN
  GIVEN( "other inputs" ){
    std::fill( wts.begin(), wts.end(), 0.0 );
    std::fill( bes.begin(), bes.end(), 0.0 );
    wts[0] = 0.8; wtn[0] = 1.0;
    int n = 0, nn = 1;
    double normEnergy = 2.5;
    std::fill( bminus.begin(), bminus.end(), 0.0 );
    bminus[0] = 1.0; bminus[1] = 2.0;
    for ( size_t i = 2; i < 10; ++i ){ bminus[i] = bminus[i-1] + bminus[i-2];}
    
    posNegTerms( n, normEnergy, bminus, wts, wtn, bes, ben, nn, -1 );
    std::vector<double> besVals {0.0, -2.5, -5.0, -7.5, -10.0, -12.5, -15.0, 
      -17.5, -20.0, -22.5, -25.0};
    std::vector<double> wtsVals {0.8, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 
      34.0, 55.0, 89.0};
    
    THEN( "values are correct" ){
      for ( size_t i = 0; i < bes.size(); ++i ){ 
        if ( i < besVals.size() ){ equal3( bes[i], besVals[i] ); }
        else { equal3( bes[i], 0.0 ); }
      }
      for ( size_t i = 0; i < wts.size(); ++i ){
        if ( i < wtsVals.size() ){ equal3( wts[i], wtsVals[i] ); }
        else { equal3( wts[i], 0.0 ); }
      }

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
      std::vector<double> besVals {0.0, -2.5, -5.0, -7.5, -10.0, -12.5, -15.0, 
        -17.5, -20.0, -22.5, -25.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
        20.0, 22.5, 25.0};
      std::vector<double> wtsVals {0.8, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 
        34.0, 55.0, 89.0, 8.0E-2, 6.0E-2, 4.0E-2, 2.0E-2, 8.0E-3, 6.0E-3, 
        4.0E-3, 2.0E-3, 8.0E-4, 6.0E-4};

      for ( size_t i = 0; i < bes.size(); ++i ){ 
        if ( i < besVals.size() ){ equal3( bes[i], besVals[i] ); }
        else { equal3( bes[i], 0.0 ); }
      }
      
      for ( size_t i = 0; i < wts.size(); ++i ){
        if ( i < wtsVals.size() ){ equal3( wts[i], wtsVals[i] ); }
        else { equal3( wts[i], 0.0 ); }
      }
      
      equal3( n, 20 );

    } // THEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "oscillator loop" ){
  GIVEN( "inputs" ){
    int maxdd = 500, numOscillators = 2, a = 0;

    double scaling = 1.0, wt = 2.0, tbart = 405.894676, temp = 200.0;

    std::vector<double> wts (maxdd, 0.0), bes(maxdd, 0.0), energyNorm(50, 0.0), 
      dbw(50, 0.0), ar(50, 0.0), dist(50,0.0);

    energyNorm[0] = 2.030778; energyNorm[1] = 2.901112;
    dist[0] = 4.55739E-3; dist[1] = 2.232634E-2;
    dbw[0] = 0.128237; dbw[1] = 0.307831;
    ar[0] = 8.213274E-2; ar[1] =  0.1368162;

    std::vector<double> alpha {0.1, 0.2, 0.4, 0.8, 1.6};
    std::vector<double> weight{0.2, 0.8};

    oscillatorLoop( alpha, dbw, ar, scaling, wts, bes, energyNorm, 
      a, maxdd, numOscillators, wt, tbart, weight, dist, temp );

    THEN( "ouput is correct" ){
      std::vector<double> wtsCorrect {0.9573911, 1.085300E-2, 6.151529E-5, 
        2.324478E-7, 1.424275E-3, 1.0594276E-6, 2.793543E-2, 3.166766E-4, 
        1.794936E-6, 4.155853E-5, 3.091273E-8, 4.075663E-4, 4.620180E-6, 
        2.618737E-8, 6.063216E-7, 3.964163E-6, 4.493784E-8, 2.891790E-8, 
        1.535389E-3, 1.740520E-5, 9.865342E-8, 2.284142E-6, 1.231187E-6, 
        1.395677E-8};
      std::vector<double> besCorrect {0.0, -2.030778, -4.061556, -6.092335, 
        2.030778, 4.061556, -2.901112, -4.931890, -6.962669, -0.8703336, 
        1.160444, -5.802224, -7.833002, -9.863781, -3.771445, -8.703336, 
        -10.73411, -11.60444, 2.901112, 0.8703336, -1.160444, 4.931890, 
        5.802224, 3.771445};

      for ( size_t i = 0; i < wts.size(); ++i ){ 
        if ( i < wtsCorrect.size() ){ equal3( wts[i], wtsCorrect[i] ); } 
        else { equal3( wts[i], 0.0 ); }
      }

      for ( size_t i = 0; i < bes.size(); ++i ){ 
        if ( i < besCorrect.size() ){ equal3( bes[i], besCorrect[i] ); } 
        else { equal3( bes[i], 0.0 ); }
      }

      equal3( wt, 3.0 );
      equal3( tbart, 407.4545311 );

    } // THEN
  } // GIVEN
} // TEST CASE
















