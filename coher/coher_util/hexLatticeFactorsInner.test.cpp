#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "hexLatticeFactorsInner.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( std::abs(b-a) < 1e-7 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-7 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}


TEST_CASE( "Function to Compute Hex Lattice Factors" ){
  double a = 1e-9, c1 = 1.5e15, c2 = 2.5e15, tsqx = 9.6e17,
    t2 = 3.5e-5, ulim = 9.6e19, c = 3.58e-8, tsq = 0, wint = 0;
  int i = 0, ifl = 1, lat = 3, nw = 60000, imax = 5;
  std::vector<double> b (60, 0.0);

  GIVEN( "l1 = 0, l2 = 0" ){
    WHEN( "only doing very few iterations (i3m is small)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 0, i3m = 3;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 98696046700994448, 0, 98696046700994448, 0,
          3.947841868E+17, 3.749690408e-8, 3.947841868E+17, 3.749690408e-8 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
    WHEN( "only doing very few iterations (i3m is small)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 0, i3m = 5;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 98696046700994448, 0, 98696046700994448, 0,
          3.947841868E+17, 3.749690408e-8, 3.947841868E+17, 3.749690408e-8,
          8.882644203E+17, 0, 8.882644203E+17, 0, 1.5791367E18, 1.4960564446E-9 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
    WHEN( "doing more iterations (i3m is moderate)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 0, i3m = 15;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 98696046700994448, 0, 98696046700994448, 0,
          3.947841868E+17, 3.749690408e-8, 3.947841868E+17, 3.749690408e-8,
          8.882644203E+17, 0, 8.882644203E+17, 0, 1.5791367E18, 1.4960564446E-9,
          2.467401167E+18, 0, 3.553057681E+18, 2.499793605E-8, 4.836106288E+18, 
	  0, 6.316546988E18, 3.67488758E-8, 7.99437978278E18, 0, 9.8696046E18, 
	  1.49987616E-8, 1.1942221650820E+19, 0, 1.42122307E19, 4.98685482E-10, 
	  1.667963189E19, 0, 1.934442515E19, 1.0713401166E-8 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
  } // GIVEN



  GIVEN( "l1 = 0, l2 = 1" ){
    WHEN( "only doing very few iterations (i3m is small)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 1, i3m = 3;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 59217626960855968., 9.488519031E-008, 
          59217626960855968., 9.488519031E-008, 1.57913673E+17, 5.71129194E-008,
          1.57913673E+17, 5.71129194E-8, 4.54001817E+17, 3.49660433E-8, 
	  4.54001817E+17, 3.49660433E-8 };

        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN

    /*
    WHEN( "only doing very few iterations (i3m is small)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 0, i3m = 5;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 98696046700994448, 0, 98696046700994448, 0,
          3.947841868E+17, 3.749690408e-8, 3.947841868E+17, 3.749690408e-8,
          8.882644203E+17, 0, 8.882644203E+17, 0, 1.5791367E18, 1.4960564446E-9 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
    
    WHEN( "doing more iterations (i3m is moderate)" ){
      THEN( "not many values are changed" ){

        int l1 = 0, l2 = 0, i3m = 15;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 98696046700994448, 0, 98696046700994448, 0,
          3.947841868E+17, 3.749690408e-8, 3.947841868E+17, 3.749690408e-8,
          8.882644203E+17, 0, 8.882644203E+17, 0, 1.5791367E18, 1.4960564446E-9,
          2.467401167E+18, 0, 3.553057681E+18, 2.499793605E-8, 4.836106288E+18, 
	  0, 6.316546988E18, 3.67488758E-8, 7.99437978278E18, 0, 9.8696046E18, 
	  1.49987616E-8, 1.1942221650820E+19, 0, 1.42122307E19, 4.98685482E-10, 
	  1.667963189E19, 0, 1.934442515E19, 1.0713401166E-8 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
    */
  } // GIVEN

  GIVEN( "l1 = 1, l2 = 1" ){
    WHEN( "only doing very few iterations (i3m is small)" ){
      THEN( "not many values are changed" ){

        int l1 = 1, l2 = 1, i3m = 3;
        hexLatticeFactorsInner( a, c1, c2, lat, nw, tsqx, b, ifl, 
            i, wint, t2, ulim, l1, l2, i3m );

        std::vector<double> bVals { 1.7765288088256790E+017, 2.1912796067443682E-007, 
		59217626960855968., 9.4885190311770216E-008, 2.7634892758356234E+017,
		0, 1.5791367366185040E+017, 5.7112919499825706E-008, 5.7243706768654573E+017,
		1.2455793604478343E-007, 4.5400181376483379E+017,3.4966043390872698E-008 };

        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
        for ( auto i = bVals.size(); i < b.size(); ++i ){ equal( b[i], 0 ); }

      } // THEN
    } // WHEN
  } // GIVEN


} // TEST CASE



