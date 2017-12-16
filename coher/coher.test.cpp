#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../catch.hpp"
#include "coher.h"

void equal( double a, double b ){
  std::cout << a << "    " << b << std::endl;

  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
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




TEST_CASE( "Function to Compute Hexagonal Lattice Factors" ){
  std::vector<double> b ( 60000, 0.0 );
  double a = 2.0e-8, tsq = 0, c1 = 1.5e15, c2 = 2.5e15, tau = 0, tsqx = 9.6e17,
    f = 0, eps = 5.0e-3, wint = 0, twopis = 39.5, t2 = 3.5e-5, ulim = 9.6e19,
    c = 3.58e-8;
  int lat = 3, nw = 60000, ifl = 1, i = 1, imax = 5;
  hexagonalLatticeFactors( a, tsq, c1, c2, lat, tau, nw, tsqx, b, ifl, f, eps,
    i, wint, twopis, t2, ulim, imax, c );
  std::vector<double> correctB { 98750002691571712., 0, 98750002691571712., 0, 3.95000010E+17, 3.18222909E-9, 3.95000010E+17, 3.18222909E-9, 8.88750024E+17, 0, 8.88750024E+17, 0, 1.58000004E+18, 8.74545088E-8, 2.46875006E+18, 6.39627061E-8, 3.55500009E+18, 8.34038989E-8, 4.83875013E+18, 9.42670649E-8, 6.32000017E+18, 1.45045490E-7, 7.99875021E+18, 1.36244484E-7, 9.87500026E+18, 2.29858595E-7, 1.19487503E+19, 2.51686065E-7, 1.42200003E+19, 3.02317783E-7, 1.66887504E+19, 3.03927633E-7, 1.93550005E+19, 4.30487131E-7, 2.22187506E+19, 4.59005787E-7, 2.52800006E+19, 5.46933463E-7, 2.85387507E+19, 5.67083424E-7, 3.19950008E+19, 6.41874124E-7, 3.56487509E+19, 7.36246592E-7, 3.95000010E+19, 8.26049002E-7, 4.35487511E+19, 7.84768579E-7, 4.77950013E+19, 7.81999073E-7, 5.22387514E+19, 7.29594243E-7, 5.68800015E+19, 7.34457572E-7, 6.17187516E+19, 7.14148462E-7, 6.67550018E+19, 7.28616470E-7, 7.19887519E+19, 7.17368653E-7, 7.74200021E+19, 7.38847206E-7, 8.30487522E+19, 7.45546030E-7, 8.88750024E+19, 7.27661629E-7, 9.48987525E+19, 1.72559207E-7, 59250000554622976. };
  //for ( auto i = 0; i < correctB.size(); ++i ){
  //  equal( b[i], correctB[i] );
 // }



} // TEST CASE




TEST_CASE( "coher" ){
  REQUIRE( true );
} // TEST CASE
