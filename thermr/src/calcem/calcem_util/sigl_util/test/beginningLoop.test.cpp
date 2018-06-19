#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/sigl_util/beginningLoop.h"

TEST_CASE( "110 120 130" ){
  std::vector<double> x(20,0.0);
  std::vector<double> y(20,0.0);
  std::vector<double> s(65,0.0);
  x[0] = 1.0; x[1] = 0.99; x[2] = -1.0;
  y[0] = 1.35700e5; y[1] = 1.35701e5; y[2] = 1.35809e5;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 };
  std::vector<double> beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
    

  // Initialize S(a,b)
  std::vector<std::vector<double>> sab(alpha.size(), 
      std::vector<double>(beta.size(),0));
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  int lasym = 0, lat = 1, iinc = 2, nlmax = 65, nl = 10, i = 3, 
      j = 0, nbin = 8;

  double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, bbm = 0.0, az = 11.9,
    tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2, tolin = 5e-2, sigmin = 1.0e-32, 
    s1bb = 1.1369180380, tol = 2.5e-2, xtol = 1.0e-5, 
    seep = 91200.0, yl = 13500, ymax = 13500, fract = 0.0, xl = -1.0, 
    eps = -1.0e-3;


  GIVEN( "inputs 1" ){

    THEN( "110-->110, 110-->120, 120-->110, 120-->130, not many iterations" ){


      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, /*az2,*/ lasym, teff, /*teff2,*/ lat, cliq, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      //ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, 
      //  bbm, az, /*az2,*/ lasym, teff, /*teff2,*/ lat, cliq, sb, sb2, iinc, eps, seep, s1bb );

      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 3 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 9 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE(-1 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 30174.6306224 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 135700.0 == Approx(yl).epsilon(1e-6) );
      //REQUIRE( 135829.6496457 == Approx(ymax).epsilon(1e-6) );

      std::vector<double> correctX = { 1.0, 0.99, -1.0, -0.005, -1.0 },
        correctY = { 135757.913, 135758.3455, 135829.6496, 135797.21918, 135809.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        if ( i < 5 ){ REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) ); }
        else        { REQUIRE( 0.0         == Approx(x[i]).epsilon(1e-6) ); }
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        if ( i < 5 ){ REQUIRE( correctY[i] == Approx(y[i]).epsilon(1e-6) ); }
        else        { REQUIRE( 0.0         == Approx(y[i]).epsilon(1e-6) ); }
      }

      REQUIRE( 271571.6756021== Approx(s[0]).epsilon(1e-6) );
      for ( size_t i = 1; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN



  GIVEN( "inputs 2" ){
    y[0] = 1.00000e5; y[1] = 1.00001e5; y[2] = 1.00009e5;

    e = 1.0e-3; ep = 1.2e-2; tev = 1.5e-1;
      tevz = 2.2e-1; 
      teff = 6.14e-0; 


    THEN( "110-->110, 110-->120, 120-->110, 120-->120, 120-->130, many iterations" ){


      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, /*az2,*/ lasym, teff, /*teff2,*/ lat, cliq, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      //ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, 
      //  bbm, az, /*az2,*/ lasym, teff, /*teff2,*/ lat, cliq, sb, sb2, iinc, eps, seep, s1bb );

      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 3 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 9 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE(-1 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 120.57844468516407 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 100000.0 == Approx(yl).epsilon(1e-6) );
      //REQUIRE( 538.71588696219601 == Approx(ymax).epsilon(1e-6) );

      std::vector<double> correctX = { 1.0, 0.99, -1.0, 0.99125, 0.990625, 
        0.9903125, 0.99015625, 0.99007813, 0.99003907, 0.99001954, 0.99000977, 
        0.99, -0.99902832, -0.99951416, -0.99975708, -0.99987854, -0.99993927, 
        -0.99996963, -0.99998481, -1.0 },
      correctY = { 461.24225912, 464.2445488, 538.7158869, 463.8735389, 
        464.0591955, 464.151910, 464.198238, 464.2213947, 464.23297, 464.238758, 
        464.2416537, 100001.0, 538.74710638, 538.7314976, 538.72369255, 
        538.7197898, 538.71783, 538.7168628, 538.71637506, 100009.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) );
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        REQUIRE( correctY[i] == Approx(y[i]).epsilon(1e-6) );
      }

      REQUIRE( 1085.2060021664768 == Approx(s[0]).epsilon(1e-6) );
      for ( size_t i = 1; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN

  GIVEN( "inputs 2" ){
    x[0] = 1.00; x[1] = 0.99; x[2] = -1.00;
    y[0] = 0.00; y[1] = 0.00; y[2] = 0.00;

    e = 1.0e-3; ep = 1.2e-2; tev = 1.5e-5;
      tevz = 2.2e-1; 
      teff = 6.14e-3; 
      xl =-1.0;
      yl = 0.0;
      i = 3;
      ymax = 1e-3;


    THEN( "110-->110, 110-->120, 120-->110, 120-->120, 120-->130, many iterations" ){


      //auto out = do_110_120_130_for_sigl(i, x, y, e, ep, tev, tevz, alpha, beta, sab, bbm, az, 
      //az2, lasym, teff, teff2, lat, cliq, sb, sb2, iinc, nl, sigmin, s, 
      //  nbin, fract, xl, j, ymax, eps, seep, yl, s1bb, tol, xtol);

      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, /*az2,*/ lasym, teff, /*teff2,*/ lat, cliq, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 0 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 8 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE( 1.0 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 0.0 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 0 == Approx(yl).epsilon(1e-6) );
      REQUIRE( 1e-3 == Approx(ymax).epsilon(1e-6) );
      

      std::vector<double> correctX = { 1.0, 0.99, 0.4925, -5.0E-3, -1.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        if ( i < correctX.size() ){
          REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) );
        }
        else {
          REQUIRE( 0.0 == Approx(x[i]).epsilon(1e-6) );
        }
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(y[i]).epsilon(1e-6) );
      }

      for ( size_t i = 0; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN

} // TEST CASE


