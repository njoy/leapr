#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "hexLatticeFactors.h"

void equal( double a, double b ){
  std::cout << a << "     " << b << std::endl;
  if (b == 0.0){ 
    REQUIRE( std::abs(b-a) < 1e-6 );
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


TEST_CASE( "Function to Compute Hex Lattice Factors" ){
  double a = 1e-9, c1 = 1.5e15, c2 = 2.5e15, tsqx = 9.6e17,
    t2 = 3.5e-5, ulim = 9.6e19, c = 3.58e-8, tsq = 0, wint = 0;
  int i = 0, ifl = 1, lat = 3, nw = 60000, imax = 5;
  std::vector<double> b (60000, 0.0);

  GIVEN( "few iterations" ){
    THEN( "outputs" ){
      hexLatticeFactors( a, tsq, c1, c2, lat, nw, tsqx, b, ifl, 
          i, wint, t2, ulim, imax, c );
      std::vector<double> b_0_100{ 98696046700994448., 0, 98696046700994448.,
        0, 3.94784186e17, 3.7496904e-8, 3.94784186e17, 3.7496904e-8,
        8.8826442e17, 0, 8.8826442e17, 0, 1.579136747e18, 1.880077432e-7,
	2.46740116e18, 2.2396046e-7, 3.55305768e18, 3.3377213359e-7,
	4.836106288e18, 6.100306e-8 };
      //std::cout << b[13] << "      " << b_0_100[13] << std::endl;
      //for( auto i = 0; i < b_0_100.size(); ++i ){ equal( b[i], b_0_100[i] ); }
      //REQUIRE( true );

    } // THEN
  } // GIVEN

  /*
  GIVEN( "inputs" ){
    a = 2e-8;
    THEN( "outputs" ){
      //hexLatticeFactors( a, tsq, c1, c2, lat, nw, tsqx, b, ifl, 
      //    i, wint, t2, ulim, imax, c );
      std::vector<double> b_0_100{ 98750002691571712., 0, 98750002691571712., 
        0, 3.950E+17, 3.748665E-8, 3.950E+17, 3.748665E-8, 8.887500E+17, 0, 
        8.887500E+17, 0, 1.580E+18, 1.324736E-6, 2.468750E+18, 1.024226E-6, 
        3.555E+18, 1.020262E-6, 4.838750E+18, 1.118530E-6, 6.320E+18, 
        1.779397E-6, 7.998750E+18, 1.491388E-6, 9.875E+18, 2.628514E-6, 
        1.194875E+19, 2.947118E-6, 1.422E+19, 3.003999E-6, 1.668875E+19, 
        3.998702E-6, 1.935500E+19, 5.188662E-6, 2.221875E+19, 5.459435E-6, 
        2.528E+19, 6.705006E-6, 2.853875E+19, 6.754011E-6, 3.199500E+19, 
        7.531615E-6, 3.564875E+19, 8.681841E-6, 3.950E+19, 9.398226E-6, 
        4.354875E+19, 9.487248E-6, 4.779500E+19, 9.381266E-6, 5.223875E+19, 
        8.342672E-6, 5.688E+19, 8.971466E-6, 6.171875E+19, 8.207045E-6, 
        6.675500E+19, 8.521837E-6, 7.198875E+19, 8.702413E-6, 7.742E+19, 
        8.318730E-6, 8.304875E+19, 9.076695E-6, 8.887500E+19, 8.668040E-6, 
        9.489875E+19, 2.048033E-6, 59250000554622976., 9.485926E-8, 
        59250000554622976., 9.485926E-8, 1.580E+17, 5.709731E-8, 1.580E+17, 
        5.709731E-8, 4.542500E+17, 3.495648E-8, 4.542500E+17, 3.495648E-8, 
        9.480E+17, 1.218753E-7, 9.480E+17, 1.218753E-7, 2.370E+17, 
        4.742963E-8, 2.370E+17, 4.742963E-8, 3.357500E+17, 3.916843E-8, 
        3.357500E+17, 3.916843E-8, 6.320E+17, 2.963580E-8, 6.320E+17, 
        2.963580E-8, 1.125750E+18, 4.664748E-7, 1.817E+18, 4.185890E-7 };

      //std::cout << b[13] << "      " << b_0_100[13] << std::endl;
      //for ( auto i = 0; i < 100; ++i ){ equal( b[i], b_0_100[i] ); }
      REQUIRE( true );

    } // THEN
  } // GIVEN
  */
} // TEST CASE

