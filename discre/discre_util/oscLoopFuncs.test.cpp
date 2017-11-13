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
  GIVEN( "inputs" ){
    int i = 0, n = 0;
    int maxdd = 500;
    double wtsn = 0.987274731;
    double besn = 0.0;
    std::vector<double> bminus (50, 0.0);
    std::vector<double> bes (maxdd, 0.0);
    std::vector<double> wts (maxdd, 0.0);
    wts[0] = 0.9872747318;
    std::vector<double> wtn (maxdd, 0.0);
    std::vector<double> ben (maxdd, 0.0);
    wtn[0] = 1.0;
    int nn = 1;
    bminus[0] = 1.119176E-2; 
    bminus[1] = 6.343540E-5; 
    bminus[2] = 2.397033E-7; 
    bminus[3] = 6.793257E-10; 
    bminus[4] = 1.540182E-12; 
    bminus[5] = 2.909946E-15; 
    bminus[6] = 4.712496E-18; 
    bminus[7] = 6.677675E-21; 
    bminus[8] = 8.410987E-24;
    double normalizedEnergy = 2.03077847;

    negativeTerms( i, n, normalizedEnergy, bminus, maxdd, wts, wtn, bes, ben, nn );
    
    equal( n, 3 );

    std::vector<double> besVals {0.0, -2.030778, -4.061556, -6.092335};
    for ( auto i = 0; i < bes.size(); ++i ){ 
      if ( i < 4 ){ equal( bes[i], besVals[i] ); }
      else { equal( bes[i], 0.0 ); }
    }

    std::vector<double> wtsVals {0.9872747, 1.119176E-2, 6.343540E-5, 2.397033E-7};
    for ( auto i = 0; i < wts.size(); ++i ){
      if ( i < 4 ){ equal( wts[i], wtsVals[i] ); }
      else { equal( wts[i], 0.0 ); }
    }

  } // GIVEN
} // TEST CASE
