#include "catch.hpp"
#include "coh/coh_util/sigcoh_util/legndr.h"


TEST_CASE( "legndr" ){
  std::vector<double> p(5), p2(5);
  double x; int np;
  GIVEN( "np is too small ( < 2 )" ){
    x = 1.5; np = 0;
    p = { 0.1, 0.2, 0.3, 0.4, 0.5 }; p2 = { 1.0, 1.5, 0.3, 0.4, 0.5 };
    legndr( x, p, np );
    for ( size_t i = 0; i < p.size(); ++i ){
      REQUIRE( p2[i] == Approx( p[i] ).epsilon(1e-6) );
    }
  } // GIVEN
  GIVEN( "np large enough ( >= 2 )" ){
    WHEN( "np == 2" ){
      x = 1.5; np = 2;
      p = { 0.1, 0.2, 0.3, 0.4, 0.5 }; p2 = { 1.0, 1.5, 2.875, 0.4, 0.5 };
      legndr( x, p, np );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p2[i] == Approx( p[i] ).epsilon(1e-6) );
      }
    } // WHEN
    WHEN( "np == 3" ){
      x = 1.5; np = 3;
      p = { 0.1, 0.2, 0.3, 0.4, 0.5 }; p2 = { 1.0, 1.5, 2.875, 6.1875, 0.5 };
      legndr( x, p, np );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p2[i] == Approx( p[i] ).epsilon(1e-6) );
      }
    } // WHEN
    WHEN( "np == 4" ){
      x = 1.5; np = 3;
      p = { 0.1, 0.2, 0.3, 0.4, 0.5 }; p2 = { 1.0, 1.5, 2.875, 6.1875, 0.5 };
      legndr( x, p, np );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p2[i] == Approx( p[i] ).epsilon(1e-6) );
      }
    } // WHEN

  } // GIVEN
} // TEST CASE
