#include <iostream>
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "coldh/coldh_util/betaLoop_util/jprimeLoop.h"


TEST_CASE( "jprime loop" ){
  std::vector<double> bex(11), rdbex(11), betan(5), sex(11);
  int j, nbx;
  double wt, be, al, x, swe, pj, y, tbart, total, out, swo;


  bex   = {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0};
  rdbex = {1.666667, 3.3333, 6.6667, 20, 5, 20, 6.6667, 3.33333, 1.66667, 0, 0};
  betan = {0.1, 0.15, 0.3, 0.6, 1.2};

  GIVEN( "jprime loop is over even values" ){

    WHEN( "molecular translations are assumed to not be free (free = false)" ){

      wt = 2.3; be = -1.2; al = 0.1; tbart = 950;

      j = 1; nbx = 10;
      x = 0.8; swe = 0.9; pj = 0.4; y = 0.3; 

      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }

       out = jPrime( j, be, x, swe, pj, bex, rdbex, sex, betan, al*wt,
         tbart, y, nbx, false, false );

       THEN( "output matches expectation" ){
         REQUIRE( out   == Approx(0.2353529421).epsilon(1e-6) );
       } // THEN

       j = 1; nbx = 10;
       x = 0.85; swe = 0.87; pj = 0.48; y = 0.35;

       sex = {5,4,3,2,1,1,1,2,2,1,5};

       out = jPrime( j, be, x, swe, pj, bex, rdbex, sex, betan, al*wt,
         tbart, y, nbx, false, false );

       THEN( "output matches expectation" ){
         REQUIRE( out   == Approx(8.87488195E-2).epsilon(1e-6) );
       } // THEN
    } // WHEN


    WHEN( "molecular translations are assumed to be free (free = true)" ){

      j = 4; be = 1.23; x = 0.1; swe = 0.2; pj = 2.4; 
         al = 2.0; wt = 0.05; y = 0.3;

       out = jPrime( j, be, x, swe, pj, bex, rdbex, sex, betan, al*wt,
         tbart, y, nbx, false, true );

       THEN( "output matches expectation" ){
         REQUIRE( out == Approx(1.9988998E-2).epsilon(1e-6) );
       } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "jprime loop for odd values" ){

    wt = 2.3; be = -1.2; al = 0.1; tbart = 950.0;
    j = 1, nbx = 10;

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      
      x = 0.8; swo = 0.9; pj = 0.4; y = 0.3;

      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = i + 1; }
      out= jPrime( j, be, x, swo, pj, bex, 
      rdbex, sex, betan, al*wt, tbart, y, nbx, true, false );

      THEN( "output matches expectation" ){
        REQUIRE( out == Approx(1.397417419).epsilon(1e-6) );
      } // THEN

      sex = { 5, 4, 3, 2, 1, 1, 1.72141595, 2.222455, 2.195246, 1.505971, 5 };
      x = 0.85; swo = 0.87; y = 0.35; pj = 0.48;
      out = jPrime( j, be, x, swo, pj, bex, rdbex, sex, betan, al*wt, 
        tbart, y, nbx, true, false );

      THEN( "output matches expectation" ){
        REQUIRE( out == Approx(8.01757644).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "molecular translations are assumed to not be free (free = false)" ){
      for ( size_t i = 0; i < sex.size(); ++i ){ sex[i] = 2 * (i + 1); }
      wt = 2.05; y = 0.35; swo = 0.47; 
      out = jPrime( j, be, x, swo, pj, bex, rdbex, sex, betan, al*wt, 
        tbart, y, nbx, true, true );

      THEN( "output matches expectation" ){
        REQUIRE( out == Approx(0.16137096895).epsilon(1e-6) );
      } // THEN
    } // WHEN

      /*
    WHEN( "debugging test" ){
      j = 1; be = -1.7615664563768665; x = 0.85293237512318631; swo =  0.67315854597803249;
      pj = 0.48411020780803327; double alphaWgt = 0.14679720469807220; tbart = 965.46997070312500;
      bex = {-1.7615664563768665,-0.88078322818843324,-0.44039161409421662,-0.22019580704710831,-0.14679720469807223, 0.14679720469807223, 0.22019580704710831, 0.44039161409421662, 0.88078322818843324,  1.7615664563768665,  0.0000000000000000};
     sex = {21.000000000000000,  16.000000000000000,  11.000000000000000,  6.0000000000000000,  1.0000000000000000,  1.0000000000000000,  4.8141700470590534,  7.0816268239446174,  6.6313307072544694,  3.6072870544081939,  0.0000000000000000};
     rdbex = { 1.1353531357048749,  2.2707062714097499,  4.5414125428194998,  13.624237628458502,  3.4060594071146242,  13.624237628458502,  4.5414125428194998,  2.2707062714097499,  1.1353531357048749,  0.0000000000000000,  0.0000000000000000};
     betan = { 0.14679720469807223, 0.22019580704710831, 0.44039161409421662, 0.88078322818843324,  1.7615664563768665};
     nbx = 10;
       y = 0.41493156681849508; 
      out = jPrime( j, be, x, swo, pj, bex, rdbex, sex, betan, alphaWgt, 
        tbart, y, nbx, true, false );

      //std::cout << out << std::endl;
      THEN( "output matches expectation" ){
     //   REQUIRE( out == Approx(0.16137096895).epsilon(1e-6) );
      } // THEN
    } // WHEN
    */
  } // GIVEN
} // TEST CASE
