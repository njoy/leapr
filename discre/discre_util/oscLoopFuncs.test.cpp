#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "oscLoopFuncs.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( std::abs(b-a) < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}


TEST_CASE( "negative n terms loop" ){
  int maxdd = 500;
  std::vector<double> bminus (50, 0.0), bes (maxdd, 0.0), wts (maxdd, 0.0),
                      wtn (maxdd, 0.0), ben (maxdd, 0.0);
 GIVEN( "inputs" ){
    int i = 0, n = 0, nn = 1;
    double wtsn = 0.987274731, besn = 0.0, normEnergy = 2.03077847;
    wts[0] = 0.9872747318, wtn[0] = 1.0;

    bminus[0] = 1.119176E-2; bminus[1] = 6.343540E-5; bminus[2] = 2.397033E-7; 
    bminus[3] = 6.79325E-10; bminus[4] = 1.54018E-12; bminus[5] = 2.909946E-15; 
    bminus[6] = 4.71249E-18; bminus[7] = 6.67767E-21; bminus[8] = 8.410987E-24;

    negativeTerms( i, n, normEnergy, bminus, wts, wtn, bes, ben, nn );
    
    THEN( "values are correct" ){
      equal( n, 3 );

      std::vector<double> besVals {0.0, -2.030778, -4.061556, -6.092335};
      for ( auto i = 0; i < bes.size(); ++i ){ 
        if ( i < 4 ){ equal( bes[i], besVals[i] ); }
        else { equal( bes[i], 0.0 ); }
      }

      std::vector<double> wtsVals {0.9872747, 1.119176E-2, 6.343540E-5, 
        2.397033E-7};
      for ( auto i = 0; i < wts.size(); ++i ){
        if ( i < 4 ){ equal( wts[i], wtsVals[i] ); }
        else { equal( wts[i], 0.0 ); }
      }
    } // THEN
  } // GIVEN
  GIVEN( "other inputs" ){
    std::fill( wts.begin(), wts.end(), 0.0 );
    std::fill( bes.begin(), bes.end(), 0.0 );
    wts[0] = 0.8; wtn[0] = 1.0;
    int i = 0, n = 0, nn = 1;
    double normEnergy = 2.5;
    std::fill( bminus.begin(), bminus.end(), 0.0 );
    bminus[0] = 1.0; bminus[1] = 2.0;
    for ( auto i = 2; i < 10; ++i ){ bminus[i] = bminus[i-1] + bminus[i-2];}
    
    negativeTerms( i, n, normEnergy, bminus, wts, wtn, bes, ben, nn );
    std::vector<double> besVals {0.0, -2.5, -5.0, -7.5, -10.0, -12.5, -15.0, 
      -17.5, -20.0, -22.5, -25.0};
    std::vector<double> wtsVals {0.8, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 
      34.0, 55.0, 89.0};
    
    THEN( "values are correct" ){
      for ( auto i = 0; i < bes.size(); ++i ){ 
        if ( i < besVals.size() ){ equal( bes[i], besVals[i] ); }
        else { equal( bes[i], 0.0 ); }
      }
      for ( auto i = 0; i < wts.size(); ++i ){
        if ( i < wtsVals.size() ){ equal( wts[i], wtsVals[i] ); }
        else { equal( wts[i], 0.0 ); }
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
    positiveTerms( i, n, normalizedEnergy, bplus, wts, wtn, bes, ben, nn ); 
    
    THEN( "values are correct" ){
      std::vector<double> besVals {0.0, -2.5, -5.0, -7.5, -10.0, -12.5, -15.0, 
        -17.5, -20.0, -22.5, -25.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
        20.0, 22.5, 25.0};
      std::vector<double> wtsVals {0.8, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 
        34.0, 55.0, 89.0, 8.0E-2, 6.0E-2, 4.0E-2, 2.0E-2, 8.0E-3, 6.0E-3, 
        4.0E-3, 2.0E-3, 8.0E-4, 6.0E-4};

      for ( auto i = 0; i < bes.size(); ++i ){ 
        if ( i < besVals.size() ){ equal( bes[i], besVals[i] ); }
        else { equal( bes[i], 0.0 ); }
      }
      
      for ( auto i = 0; i < wts.size(); ++i ){
        if ( i < wtsVals.size() ){ equal( wts[i], wtsVals[i] ); }
        else { equal( wts[i], 0.0 ); }
      }
      
      equal( n, 20 );

    } // THEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "oscillator loop" ){
  GIVEN( "inputs" ){
    std::vector<double> wtsCorrect {0.9573911, 1.085300E-2, 6.151529E-5, 
      2.324478E-7, 1.424275E-3, 1.059427E-6, 2.793543E-2, 3.166766E-4, 
      1.794936E-6, 4.155853E-5, 3.091273E-8, 4.075663E-4, 4.620180E-6, 
      2.618737E-8, 6.063216E-7, 3.964163E-6, 4.493784E-8, 2.891790E-8, 
      1.535389E-3, 1.740520E-5, 9.865342E-8, 2.284142E-6, 1.231187E-6, 
      1.395677E-8};

    int maxdd = 500;
    std::vector<double> wts (maxdd, 0.0), bes(maxdd, 0.0), energyNorm(50, 0.0), 
      dbw(50, 0.0), wtn(maxdd, 0.0), ben(maxdd, 0.0), ar(50, 0.0);
    energyNorm[0] = 2.030778; energyNorm[1] = 2.901112;
    dbw[0] = 0.128237; dbw[1] = 0.307831;
    wtn[0] = 1.0;
    double scaling = 1.0; int a = 0;
    ar[0] = 8.213274E-2; ar[1] =  0.1368162;
    std::vector<double> alpha {0.1, 0.2, 0.4, 0.8, 1.6};
    oscillatorLoop( alpha,  











  } // GIVEN
} // TEST CASE
















