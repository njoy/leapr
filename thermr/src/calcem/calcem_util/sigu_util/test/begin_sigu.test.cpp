#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/sigu_util/begin_sigu.h"


TEST_CASE( "begin sigu (113,116)" ){
  GIVEN( "inputs" ){
  int jbeta = -7, lat = 1, iinc = 2, 
      lasym = 0;

  std::cout << std::setprecision(10) ;
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



    do_113_116( jbeta, lat, x, y, e, tev, tevz, root1, u, alpha, beta, 
        sab, az, lasym, az2, teff, teff2, cliq, sb, sb2, iinc);

    REQUIRE( 2.3e-2 == Approx(x[0]).epsilon(1e-6) );
    REQUIRE( 0.0    == Approx(y[0]).epsilon(1e-6) ); 

    for ( size_t i = 1; i < x.size(); ++i ){ 
      REQUIRE( 0 == Approx(x[i]).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(y[i]).epsilon(1e-6) ); 
    }
 
    REQUIRE( jbeta == 1 );

  } // GIVEN


  GIVEN( "inputs 2" ){
    int jbeta = -7, lat = 1, iinc = 2, 
        lasym = 0;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-3, tev = 1.5e-5, az = 11.9,
      tevz = 2.2e-1, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2, u = 0.05, root1 = 2.9e-2;

    std::vector<double> alpha(40), beta(80), x(20,0.0), y(20,0.0);
    for ( size_t j = 0; j < alpha.size(); ++j ){ alpha[j] = 0.1*(j+1); }
    for ( size_t j = 0; j < beta.size();  ++j ){ beta[j]  = 0.5*(j+1); }


    std::vector<std::vector<double>> sab(alpha.size(), 
        std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 5.0 + i*0.03 - j*0.004;
      } 
    } 


    do_113_116( jbeta, lat, x, y, e, tev, tevz, root1, u, alpha, beta, 
        sab, az, lasym, az2, teff, teff2, cliq, sb, sb2, iinc);

    REQUIRE( 0.111 == Approx(x[0]).epsilon(1e-6) );
    REQUIRE( 0.0   == Approx(y[0]).epsilon(1e-6) ); 

    for ( size_t i = 1; i < x.size(); ++i ){ 
      REQUIRE( 0 == Approx(x[i]).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(y[i]).epsilon(1e-6) ); 
    }
 
    REQUIRE( 1 == jbeta );

  } // GIVEN



  GIVEN( "inputs 3" ){
    int jbeta = -80, lat = 1, iinc = 2, 
        lasym = 0;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-6, tev = 1.5e-4, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2, u = 0.1, root1 = 0.00092700235;

    std::vector<double> alpha(40),beta(80),x(20,0.0), y(20,0.0);
    for ( int i = 0; i < 40;     ++i ){ alpha[i] = 0.1*i + i%6 + 0.001; }
    for ( int i = 0; i < 80;     ++i ){ beta[i]  = 0.2*i + i%4 + 0.025; }

    std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.2*i + 0.4*j + (i+j)%5;
      } 
    }

    do_113_116( jbeta, lat, x, y, e, tev, tevz, root1, u, alpha, beta, 
        sab, az, lasym, az2, teff, teff2, cliq, sb, sb2, iinc);

    REQUIRE( 6.5e-6 == Approx(x[0]).epsilon(1e-6) );
    REQUIRE( 46226.360425 == Approx(y[0]).epsilon(1e-6) ); 

    for ( size_t i = 1; i < x.size(); ++i ){ 
      REQUIRE( 0 == Approx(x[i]).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(y[i]).epsilon(1e-6) ); 
    }
 
    REQUIRE( 1 == jbeta );

  } // GIVEN



} // TEST CASE

//    std::cout << x[0] << std::endl;
//    std::cout << y[0] << std::endl;


