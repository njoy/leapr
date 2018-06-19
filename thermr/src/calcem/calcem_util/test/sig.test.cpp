#include "catch.hpp"
#include "calcem/calcem_util/sig.h"


TEST_CASE( "sig" ){
  int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc;

  double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, az = 11.9,
    tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2;

  double sigVal1, sigVal2;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
    beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i][j] = 0.01*((j+1) + 0.1*(i+1));

    } 
  } 

  GIVEN( "invalid input requested (iinc != 1,2)" ){
    iinc = 0;
    THEN( "throw exception" ){
      REQUIRE_THROWS( sig( e, ep, u, tev, alpha, beta, sab, 
        az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc ) );
    } // THEN
  } // GIVEN

  GIVEN( "free gas option selected (iinc = 1)" ){
    iinc = 1;

    WHEN( "final neutron energy E' is zero (ep = 0)" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    ep = 1.2e-6;

    WHEN( "temperature is extremely small" ){
      tev = 1.5e-8;
      sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
        az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
      THEN( "alpha, beta --> big, so exponential in Eq. 229 dies" ){
        REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "nonextreme values are provided" ){
      sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
        az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
      THEN( "nonzero sig output is provided" ){
        REQUIRE( 1384.063418 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN ("alpha, beta values outside of provided range" ){
    iinc = 2; tev = 1.5e-1; tevz = 2.2e-4;
    WHEN( "this is because of a > alpha(nalpha)" ){
      e = 1.0e-2, ep = 1.2e-2;
      sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
        az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );

      e = 1.0e-3, ep = 1.2e-3; alpha = { 0.1, 0.2, 0.3, 0.5, 0.8 };
      sigVal2 = sig( e, ep, u, tev, alpha, beta, sab, 
        az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );

      THEN( "SCT Approximation is used" ){
        REQUIRE( 55.838635 == Approx( sigVal1 ).epsilon(1e-6) );
        REQUIRE( 179.2164755 == Approx( sigVal2 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "this is because of being out of range with beta" ){
      e = 1.0e-5;

      AND_WHEN( "lasym != 1" ){ 
        e = 1.0e-4, ep = 5.2e-3, tev = 1.5e-1;
        lasym = 0;
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 795.18855771477308 == Approx( sigVal1 ).epsilon(1e-6) );
      } // AND WHEN

      AND_WHEN( "lasym == 1" ){
        AND_WHEN( "bbm is smaller than min value in beta vector" ){ 
          ep = 1.2e-6, tev = 1.5e-8;
          alpha = { 0.1, 0.2, 0.3, 0.5, 0.8 };
          lasym = 1;
          sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
            az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
          REQUIRE( 883.34692414 == Approx( sigVal1 ).epsilon(1e-6) );
        } // AND WHEN
        AND_WHEN( "bbm is greater than max value in beta vector" ){ 
            ep = 1.2e-1;
            alpha = { 10.1, 20.2, 30.3, 40.5, 50.8 };

          AND_WHEN( "-arg < sabflg" ){ 
             tev = 1.5e-8;
            sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
              az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
          } // AND WHEN

          AND_WHEN( "-arg > sabflg" ){ 
            sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
              az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );

            e = 1.0e-4;
            sigVal2 = sig( e, ep, u, tev, alpha, beta, sab, 
              az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
            THEN( "SCT Approximation is used" ){
              REQUIRE( 12.9236152 == Approx( sigVal1 ).epsilon(1e-6) );
              REQUIRE( 5.011738726 == Approx( sigVal2 ).epsilon(1e-6) );
            } // THEN
          } // AND WHEN
        } // AND WHEN
      } // AND WHEN
    } // WHEN
  } // GIVEN

  iinc = 2;
  e = 1.0e-6, ep = 1.2e-4, u = -1.0, tev = 1.5e-1, tevz = 2.2e-4;

  GIVEN( "150" ){
    WHEN( "final neutron energy E' is zero (ep = 0)" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 201.87960468 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    tev = 1.5e-4;
    WHEN( "150 to 155 to 160" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 135829.64964 == Approx( sigVal1 ).epsilon(1e-6) );

        cliq = 1.0;
        sigVal2 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 135829.64964 == Approx( sigVal2 ).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "alpha, beta values are smaller" ){
      for ( auto& entry : alpha ){ entry *= 0.01; }
      for ( size_t i = 0; i < beta.size()-1; ++i  ){ beta[i] *= 0.1; }
      THEN( "search for ia and ib in alpha, beta vectors" ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 187405.9625716 == Approx( sigVal1 ).epsilon(1e-6) );

      } // THEN
    } // WHEN

  } // GIVEN
  GIVEN( "160" ){
    tev = 1.5e-4;
    WHEN( "straight to 160" ){
      cliq = 1.0;
      ep = 3.1e-5;
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 248176.610043 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( "other stuff that's failing" ){
    WHEN( "straight to 160" ){
      cliq = 1.0; e = 1e-6; tev = 1.5e-4; ep = 2.3e-5; tevz = 2.2e-4; 
      teff = 6.14e-2; u = 0.1;
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 298101.75219413883 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "lat does not equal 1" ){
    WHEN( "" ){

      e = 1e-6; ep = 1e-3; u = 0.1;
      tev = 1.5e-4; az = 11.9; tevz = 2.2e-4;
      lasym = 0; az2 = 0; teff2 = 0;
      lat = 0; cliq = 1; sb = 5.53;
      sb2 = 0; teff = 6.14e-2; iinc = 2;
      THEN( "" ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 13.39089778 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "extra" ){
    WHEN( "" ){
    std::vector<double> alpha { 0.25203, 0.50406, 0.75609, 1.00812, 1.26015, 
      1.51218, 1.76421, 2.01624, 2.27331, 2.53552, 2.80297, 3.07577, 3.35401, 
      3.63790, 3.92733, 4.22271, 4.52383, 4.83111, 5.14443, 5.46411, 5.79013, 
      6.12261, 6.46185, 6.80783, 7.16077, 7.52067, 7.88783, 8.26234, 8.64432, 
      9.03396, 9.43136, 9.83673, 10.2506, 10.6719, 11.1024, 11.5409, 11.9886, 
      12.4452, 12.911, 13.3858 }, 
      beta { 0.0, 0.100812, 0.201624, 0.302436, 0.403248, 0.50406, 0.604872, 
      0.705684, 0.806496, 0.907307, 1.00812, 1.10893, 1.20974, 1.31055, 
      1.41137, 1.51218, 1.61299, 1.71380, 1.81461, 1.91543, 2.01624, 
      2.11705, 2.21786, 2.31867, 2.41949, 2.5203, 2.62111, 2.72192, 
      2.82273, 2.92354, 3.02436, 3.12517, 3.22598, 3.32679, 3.42760, 
      3.52842, 3.62923, 3.73004, 3.83085, 3.93167, 4.03248, 4.13329, 
      4.24378, 4.36485, 4.49762, 4.64309, 4.80248, 4.97719, 5.16873, 
      5.37862, 5.60867, 5.8608, 6.13713, 6.43997, 6.77184, 7.13567, 
      7.53438, 7.97140, 8.45036, 8.97529, 9.55052, 10.181, 10.8726, 
      11.6297, 12.4593, 13.3697, 14.3667, 15.4595, 16.6571, 17.9697, 
      19.4093, 20.9860, 22.7139, 24.6082, 26.6849, 28.9602, 31.4533, 
      34.1873, 37.1825, 40.4659 };
    std::vector<std::vector<double>> sab ( alpha.size(), std::vector<double> (beta.size()));
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.5*i + 0.1*j;
      }
    }


      
      e = 1e-5; ep = 9.999999e-06; u = 0.99;
      tev = 0.0255074596; az = 11.9; tevz = 0.0253;
      lasym = 0; az2 = 0; teff2 = 0;
      lat = 1; cliq = 0; sb = 5.53486;
      sb2 = 0; teff = 0.0614755628515; iinc = 2;
      THEN( "" ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 1.3856176085894635E-2 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "extra" ){
    WHEN( "" ){
    std::vector<double> alpha(40),beta(80);
    for ( int i = 0; i < 40; ++i ){ alpha[i] = 0.1*i + 0.001; }
    for ( int i = 0; i < 80; ++i ){ beta[i]  = 0.2*i + 0.025; }

    std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.2*i + 0.4*j + (i+j)%5;
      } 
    } 

      e = 1e-6; ep = 0.0005345; u = 0.1;
      tev = 0.00015; az = 11.9; tevz = 0.00022;
      lasym = 0; az2 = 0; teff2 = 0;
      lat = 1; cliq = 0; sb = 5.53;
      sb2 = 0; teff = 6.14e-1; iinc = 2;
      THEN( "" ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 665890150.16903055 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "extra" ){
    WHEN( "" ){
    std::vector<double> alpha(40),beta(80);
    for ( int i = 0; i < 40; ++i ){ alpha[i] = 0.1*i + 0.001; }
    for ( int i = 0; i < 80; ++i ){ beta[i]  = 0.2*i + 0.025; }

    std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.2*i + 0.4*j + (i+j)%5;
      } 
    } 

      e = 1e-6; ep = 6.5e-6; u = 0.1;
      tev = 0.00015; az = 11.9; tevz = 0.00022;
      lasym = 0; az2 = 0; teff2 = 0;
      lat = 1; cliq = 0; sb = 5.53;
      sb2 = 0; teff = 6.14e-1; iinc = 2;
      THEN( "" ){
        sigVal1 = sig( e, ep, u, tev, alpha, beta, sab, 
          az, tevz, lasym, /*az2, teff2,*/lat, cliq, sb, sb2, teff, iinc );
        //REQUIRE( 665890150.16903055 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

} // TEST CASE

