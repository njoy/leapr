#include "catch.hpp"
#include "calcem/calcem_util/sigu_util/do150.h"


TEST_CASE( "do 150" ){
  GIVEN( "inputs" ){
  int jbeta = -7, lat = 1, iinc = 2, 
      lasym = 0;

  std::cout << std::setprecision(15) ;
  double e = 1.0e-3, tev = 1.5e-5, az = 11.9,
    tevz = 2.2e-1, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2, tolin = 5e-2, u, root1 = 2.9e-2;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
    beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0), y(20,0.0);

  std::vector<std::vector<double>> sab(alpha.size(), 
      std::vector<double>(beta.size(),0));
  for ( int i = 0; i < alpha.size(); ++i ){
    for ( int j = 0; j < beta.size(); ++j ){
      sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 


  do150( i, x, y, xm, ym, yt, test, tolmin, e, u, tev, alpha, beta, sab, tevz, lasym, az, az2, teff2, lat, cliq, sb, sb2, teff, tol, iinc);

 

  } // GIVEN
} // TEST CASE


