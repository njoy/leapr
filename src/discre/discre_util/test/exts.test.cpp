#include "catch.hpp"
#include "discre/discre_util/exts.h"

void equal1( double a, double b ){
	std::cout << a << "       " << b << std::endl;
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal1_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal1( a[i], b[i] );
  }
}



TEST_CASE( "exts" ){
  std::vector<double> beta(5), sexpb(5), exb(11), sex(11), correctSex(11);

  GIVEN( "inputs" ){

    WHEN( "First beta value is really small (<=1e-9)" ){
      beta = { 1.0e-9, 0.15, 0.30, 0.60, 1.20 };
      sexpb = { 53, 9E-2, 2E-2, 6E-3, 1E-3 };
      exb = { 1.0, 0.9, 0.8, 0.7, 0.5 };

      sex = exts( sexpb, exb, beta );
      correctSex = {  };


      equal1_vec( sex, correctSex );
      /*
      
      sex = exts( sexpb, exb, beta );
      equal1_vec( sex, correctSex );
      */

    } // WHEN
    /*
    WHEN( "First beta value is not that small (>1e-9)" ){
      beta  = { 0.1, 0.15, 0.30, 0.60, 1.20 };
      sexpb = {0.1992005, 9.052783E-2, 2.393197E-2, 6.549547E-3, 1.876351E-3};
      exb = {0.9512294, 0.9277434, 0.8607079, 0.7408182, 0.5488116, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0};
      std::vector<double> correctSex =  {1.876351E-3, 6.549547E-3, 2.393197E-2, 
        9.052783E-2, 0.1992005, 0.1992005, 7.791802E-2, 1.772923E-2, 
        3.594467E-3, 5.651460E-4, 0.0};
      auto sex = exts( sexpb, exb, beta );
      equal1_vec( sex, correctSex );

      beta  = { 1.1e-9, 0.15, 0.30, 0.60, 1.20 };
      sex = exts( sexpb, exb, beta );
      equal1_vec( sex, correctSex );
    } // WHEN 
    */
  } // GIVEN
} // TEST CASE
