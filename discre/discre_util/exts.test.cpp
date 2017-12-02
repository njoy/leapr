#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "exts.h"

void equal( double a, double b ){
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



TEST_CASE( "exts" ){
  GIVEN( "inputs" ){
    std::vector<double> beta { 0.0, 0.15, 0.30, 0.60, 1.20 };
    std::vector<double> sexpb {53.06972, 9.052783E-2, 2.393197E-2, 
      6.549547E-3, 1.876351E-3};
    std::vector<double> exb {1.0, 0.9277434, 0.8607079, 0.7408182, 0.5488116,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    WHEN( "First beta value is really small (<=1e-9)" ){
      auto sex = exts( sexpb, exb, beta );
      std::vector<double> correctSex {1.876351E-3, 6.549547E-3, 2.393197E-2, 
        9.052783E-2, 53.06972, 7.791802E-2, 1.772923E-2, 3.594467E-3, 
        5.651460E-4, 0.0, 0.0};
      equal_vec( sex, correctSex );
      
      beta = { 1.0e-9, 0.15, 0.30, 0.60, 1.20 };
      sex = exts( sexpb, exb, beta );
      equal_vec( sex, correctSex );

    } // WHEN
    WHEN( "First beta value is not that small (>1e-9)" ){
      beta  = { 0.1, 0.15, 0.30, 0.60, 1.20 };
      sexpb = {0.1992005, 9.052783E-2, 2.393197E-2, 6.549547E-3, 1.876351E-3};
      exb = {0.9512294, 0.9277434, 0.8607079, 0.7408182, 0.5488116, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0};
      std::vector<double> correctSex =  {1.876351E-3, 6.549547E-3, 2.393197E-2, 
        9.052783E-2, 0.1992005, 0.1992005, 7.791802E-2, 1.772923E-2, 
        3.594467E-3, 5.651460E-4, 0.0};
      auto sex = exts( sexpb, exb, beta );
      equal_vec( sex, correctSex );

      beta  = { 1.1e-9, 0.15, 0.30, 0.60, 1.20 };
      sex = exts( sexpb, exb, beta );
      equal_vec( sex, correctSex );
    } // WHEN 
  } // GIVEN
} // TEST CASE
