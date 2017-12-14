#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../catch.hpp"
#include "coher.h"


TEST_CASE( "Function to Compute Hexagonal Lattice Factors" ){
  std::vector<double> b ( 60000, 0.0 );
  double a = 2.0e-8, tsq = 0, c1 = 1.5e15, c2 = 2.5e15, tau = 0, tsqx = 9.6e17,
    f = 0, eps = 5.0e-3, wint = 0, twopis = 39.5, t2 = 3.5e-5, ulim = 9.6e19,
    c = 3.58e-8;
  int lat = 3, nw = 60000, ifl = 1, i = 1, imax = 5;

  hexagonalLatticeFactors( a, tsq, c1, c2, lat, tau, nw, tsqx, b, ifl, f, eps,
    i, wint, twopis, t2, ulim, imax, c );
  std::cout << b[0] << std::endl;



} // TEST CASE




TEST_CASE( "coher" ){
  REQUIRE( true );
} // TEST CASE
