#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/end.h"

void equal1( double a, double b ){
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

void equal1Vec( const std::vector<double>& a, const std::vector<double>& b ){
  REQUIRE( a.size() == b.size() );
  for ( size_t i = 0; i < a.size(); ++i ){
    equal1( a[i], b[i] );
  }
}


TEST_CASE( "sortLatticeFactors" ){
  std::vector<double> b(20), correctSort(20);
  GIVEN( "inputs" ){
    int ifl = 1, k = 4, imax = 12;
    double ulim = 9e19;
    b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
    correctSort = {13,14,15,21,19,12,41,8,9e19,8,21,81,53,37,39,93,42,45,57,3};
    sortLatticeFactors( ifl, b, k, ulim, imax );
    equal1Vec( b, correctSort );
    
    b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
    correctSort = {13,14,15,21,41,8,19,12,9e19,12,21,81,53,37,39,93,42,45,57,3};
    imax = 2; k = 4;
    sortLatticeFactors( ifl, b, k, ulim, imax );
    equal1Vec( b, correctSort );
    
  } // GIVEN
} // TEST CASE

TEST_CASE( "end" ){
  std::vector<double> b(20), correctSort(20);
  GIVEN( "inputs" ){
    {
      int ifl = 1, k = 4, imax = 2, nbe, nw;
      double toler = 1e-6, recon = 5e-20, scon = 12000, ulim = 9e19;
      b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
      correctSort = {6.5E-19,660000,4.5,144000,41,8,19,12,9e19,12,21,81,53,37,39,93,42,45,57,3};
      nbe = end( ifl, b, k, recon, toler, scon, ulim, imax );
      equal1Vec( b, correctSort );
      nw = 2*k;
      REQUIRE( nbe == 2  );
      REQUIRE( nw  == 10 );
    }


    {
      int ifl = 1, k = 6, imax = 12, nbe, nw;
      double toler = 1e-6, recon = 1e-20, scon = 123000, ulim = 9e19;
      b = {19,12,15,21,41,8,13,14,17,0,0,0,0,0,0,0,0,0,0,0};
      correctSort = {0.0, 6765000.0, 0.9, 984000.0, 15, 21, 17, 0, 19, 12, 41, 8, 9E19, 8, 0, 0, 0, 0, 0, 0};
      nbe = end( ifl, b, k, recon, toler, scon, ulim, imax );
      equal1Vec( b, correctSort );
      nw = 2*k;
      REQUIRE( nbe == 2  );
      REQUIRE( nw  == 14 );
    }

    {
      int ifl = 1, k = 6, imax = 12, nbe, nw;
      double toler = 1e-6, recon = 1e-20, scon = 123000, ulim = 9e19;
      std::vector<double> b(200);
      for (size_t i = 0; i < b.size(); ++i){
        if (i%2 == 0){ b[i]   = i*1e-3 + 0.2;}
        if (i%2 == 1){ b[i] = i*1e3  + 200;}
      }
      std::vector<double> correctSort = {2.E-21, 4.5756E9, 9.0E-1, 1.3776E9, 
      2.04E-1, 5.2E3, 2.06E-1, 7.2E3, 2.08E-1, 9.2E3, 2.1E-1, 1.12E4, 9.E19, 
      1.12E4, 2.14E-1, 1.52E4, 2.16E-1, 1.72E4, 2.18E-1, 1.92E4, 2.2E-1, 2.12E4, 
      2.22E-1, 2.32E4, 2.24E-1, 2.52E4, 2.26E-1, 2.72E4, 2.28E-1, 2.92E4, 2.3E-1, 
      3.12E4, 2.32E-1, 3.32E4, 2.34E-1, 3.52E4, 2.36E-1, 3.72E4, 2.38E-1, 3.92E4, 
      2.4E-1, 4.12E4, 2.42E-1, 4.32E4, 2.44E-1, 4.52E4, 2.46E-1, 4.72E4, 2.48E-1, 
      4.92E4, 2.5E-1, 5.12E4, 2.52E-1, 5.32E4, 2.54E-1, 5.52E4, 2.56E-1, 5.72E4, 
      2.58E-1, 5.92E4, 2.6E-1, 6.12E4, 2.62E-1, 6.32E4, 2.64E-1, 6.52E4, 2.66E-1,
      6.72E4, 2.68E-1, 6.92E4, 2.7E-1, 7.12E4, 2.72E-1, 7.32E4, 2.74E-1, 7.52E4, 
      2.76E-1, 7.72E4, 2.78E-1, 7.92E4, 2.8E-1, 8.12E4, 2.82E-1, 8.32E4, 2.84E-1, 
      8.52E4, 2.86E-1, 8.72E4, 2.88E-1, 8.92E4, 2.9E-1, 9.12E4, 2.92E-1, 9.32E4, 
      2.94E-1, 9.52E4, 2.96E-1, 9.72E4, 2.98E-1, 9.92E4, 3.E-1, 1.012E5, 3.02E-1, 
      1.032E5, 3.04E-1, 1.052E5, 3.06E-1, 1.072E5, 3.08E-1, 1.092E5, 3.1E-1, 
      1.112E5, 3.12E-1, 1.132E5, 3.14E-1, 1.152E5, 3.16E-1, 1.172E5, 3.18E-1, 
      1.192E5, 3.2E-1, 1.212E5, 3.22E-1, 1.232E5, 3.24E-1, 1.252E5, 3.26E-1, 
      1.272E5, 3.28E-1, 1.292E5, 3.3E-1, 1.312E5, 3.32E-1, 1.332E5, 3.34E-1, 
      1.352E5, 3.36E-1, 1.372E5, 3.38E-1, 1.392E5, 3.4E-1, 1.412E5, 3.42E-1, 
      1.432E5, 3.44E-1, 1.452E5, 3.46E-1, 1.472E5, 3.48E-1, 1.492E5, 3.5E-1, 
      1.512E5, 3.52E-1, 1.532E5, 3.54E-1, 1.552E5, 3.56E-1, 1.572E5, 3.58E-1, 
      1.592E5, 3.6E-1, 1.612E5, 3.62E-1, 1.632E5, 3.64E-1, 1.652E5, 3.66E-1, 
      1.672E5, 3.68E-1, 1.692E5, 3.7E-1, 1.712E5, 3.72E-1, 1.732E5, 3.74E-1, 
      1.752E5, 3.76E-1, 1.772E5, 3.78E-1, 1.792E5, 3.8E-1, 1.812E5, 3.82E-1, 
      1.832E5, 3.84E-1, 1.852E5, 3.86E-1, 1.872E5, 3.88E-1, 1.892E5, 3.9E-1, 
      1.912E5, 3.92E-1, 1.932E5, 3.94E-1, 1.952E5, 3.96E-1, 1.972E5, 3.98E-1, 
      1.992E5};
      nbe = end( ifl, b, k, recon, toler, scon, ulim, imax );
      equal1Vec( b, correctSort );
      nw = 2*k;
      REQUIRE( nbe == 2  );
      REQUIRE( nw  == 14 );
    }




  } // GIVEN
} // TEST CASE 
    
    









