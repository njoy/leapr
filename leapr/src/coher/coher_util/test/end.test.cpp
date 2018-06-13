
#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/end.h"

void equal1( double a, double b ){
  if( b == 0 ){ REQUIRE( std::abs(a-b) < 1e-6 ); }
  if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

void equal1Vec( const std::vector<double>& a, const std::vector<double>& b ){
  REQUIRE( a.size() == b.size() );
  for ( auto i = 0; i < a.size(); ++i ){
    equal1( a[i], b[i] );
  }
}


TEST_CASE( "sortLatticeFactors" ){
  std::vector<double> b(20), correctSort(20);
  GIVEN( "inputs" ){
    int ifl = 1, k = 4, nw = 20, imax = 12;
    double ulim = 9e19;
    b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
    correctSort = {13,14,15,21,19,12,41,8,9e19,8,21,81,53,37,39,93,42,45,57,3};
    sortLatticeFactors( ifl, b, k, nw, ulim, imax );
    equal1Vec( b, correctSort );
    
    b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
    correctSort = {13,14,15,21,41,8,19,12,9e19,12,21,81,53,37,39,93,42,45,57,3};
    imax = 2; k = 4;
    sortLatticeFactors( ifl, b, k, nw, ulim, imax );
    equal1Vec( b, correctSort );
    
  } // GIVEN
} // TEST CASE


TEST_CASE( "end" ){
  std::vector<double> b(20), correctSort(20);
  GIVEN( "inputs" ){
    int ifl = 1, k = 4, nw = 20, imax = 2, maxb = 20;
    double toler = 1e-6, recon = 5e-20, scon = 12000, ulim = 9e19;
    b = {19,12,15,21,41,8,13,14,17,74,21,81,53,37,39,93,42,45,57,3};
    correctSort = {6.5E-19,660000,4.5,144000,41,8,19,12,9e19,12,21,81,53,37,39,93,42,45,57,3};

    end( ifl, b, k, recon, maxb, toler, scon, nw, ulim, imax );
    equal1Vec( b, correctSort );
    

  } // GIVEN
} // TEST CASE 
    
    









