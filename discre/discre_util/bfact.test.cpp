#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "bfact.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}



TEST_CASE( "bfact" ){
  GIVEN( "inputs" ){
    std::vector<double> bplus (50, 0.0), bminus (50, 0.0);
    double x   = 8.213274e-3, dwc = 1.282379e-2, beta_i = 2.030778478;

    double bzero = bfact(x, dwc, beta_i, bplus, bminus );

    std::vector<double> bminusStart {1.119176E-2, 6.343540E-5, 2.397033E-7, 
      6.793257E-10, 1.540182E-12, 2.909946E-15, 4.712496E-18, 6.677675E-21, 
      8.410987E-24};
    std::vector<double> bplusStart {1.468732E-3, 1.092496E-6, 5.417594E-10, 
      2.014904E-13, 5.995048E-17, 1.486447E-20, 3.159074E-24, 5.874598E-28};
    for ( auto i = 0; i < bminusStart.size(); ++i ){
      equal( bminus[i], bminusStart[i] );
   }
    for ( auto i = bminusStart.size(); i < bminus.size(); ++i ){
      equal( bminus[i], 0.0 );
    }
    for ( auto i = 0; i < bplusStart.size(); ++i ){
      equal( bplus[i], bplusStart[i] );
   }
    for ( auto i = bplusStart.size(); i < bplus.size(); ++i ){
      equal( bplus[i], 0.0 );
    }
    equal( bzero, 0.98727473 );



    x = 0.01, dwc = 0.2, beta_i = 3.0;
    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    std::fill (bplus.begin(),bplus.end(),0);
    std::fill (bminus.begin(),bminus.end(),0);

    bzero = bfact(x, dwc, beta_i, bplus, bminus );
    
    bplusStart = {9.134290E-4, 5.095317E-7, 1.894861E-10, 5.285001E-14, 
      1.179242E-17, 2.192702E-21, 3.494699E-25, 4.873578E-29};
    bminusStart = {1.834671E-2, 2.055597E-4, 1.535421E-6, 8.601593E-9, 
      3.854963E-11, 1.439728E-13, 4.608864E-16, 1.290968E-18, 3.214286E-21, 
      7.202715E-24};
    for ( auto i = 0; i < bminusStart.size(); ++i ){
      equal( bminus[i], bminusStart[i] );
    }
    for ( auto i = bminusStart.size(); i < bminus.size(); ++i ){
      equal( bminus[i], 0.0 );
    }
    for ( auto i = 0; i < bplusStart.size(); ++i ){
      equal( bplus[i], bplusStart[i] );
    }
    for ( auto i = bplusStart.size(); i < bplus.size(); ++i ){
      equal( bplus[i], 0.0 );
    }
    equal( bzero, 0.818751219 );


 } // GIVEN
} // TEST CASE
