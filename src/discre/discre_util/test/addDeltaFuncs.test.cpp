#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "discre/discre_util/addDeltaFuncs.h"


void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( size_t i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}

TEST_CASE( "helper function to add delta functions to scattering law" ){
  double twt = -0.1, dwf = 0.9;
  int n = 17;

  std::vector<double> bes ( 500, 0.0 ), wts ( 500, 0.0 );

  std::vector<double> betan { 0.1, 0.15, 0.30, 0.60, 1.20 }, 
                      sexpb { 0.1, 0.2,  0.3,  0.4,  0.5  };

    double val = 9; int sign = -1;
    for ( auto i = 1; i < 20; ++i ){ 
      bes[i] = sign * val;
      sign *= -1; 
      val -= 1;
    }
    int j = 0;
    for ( auto i = 0.01; i > 1.0e-8; i = i / 10 ){
      wts[j] = 4.0*i;
      wts[j+1] = 3.0*i;
      wts[j+2] = 2.0*i;
      wts[j+3] = 1.0*i;
      j = j + 4;
    }

  GIVEN( "all neg bes vals have larger mag. than the 2nd largeset beta val" ){
  //  addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
    THEN( "no change to sexpb" ){
      equal_vec( sexpb, {0.1,0.2,0.3,0.4,0.5} );
    } // THEN
  } // GIVEN

  GIVEN( "some neg bes vals have smaller mag. than the 2nd largeset beta val" ){
    sexpb = { 0.1, 0.2,  0.3,  0.4,  0.5  };
    betan[3] = 2.3; betan[4] = 2.5;
    addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
    THEN( "sexpb is correctly amended" ){
      equal_vec( sexpb, {0.1,0.2,0.3002511747,0.4,0.5} );
    } // THEN
  } // GIVEN

  GIVEN( "Large difference between beta values and bes values" ){
    sexpb = { 0.1, 0.2,  0.3,  0.4,  0.5  };
    betan = { 1001, 1002, 1003, 1004, 1005 };
    addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
    THEN( "reslting index jj is so small you can't access the (jj-2)th beta" ){
      equal_vec( sexpb, {0.1000269475,0.2,0.3,0.4,0.5} );
    } // THEN
    WHEN( "twt value is changed from one negative val to another" ){
      sexpb = { 0.1, 0.2,  0.3,  0.4,  0.5  };
      twt = -15.0;
      THEN( "no change in output" ){
        addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
        equal_vec( sexpb, {0.1000269475,0.2,0.3,0.4,0.5} );
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "A small (< 1.0e-10) vallue for dwt" ){
    sexpb = { 0.1, 0.2,  0.3,  0.4,  0.5  };
    betan = { 1001, 1002, 1003, 1004, 1005 };
    dwf = 9.9e-11;
    addDeltaFuncs( twt, dwf, bes, betan, wts, sexpb, n );
    THEN( "no change to sexpbc" ){
      equal_vec( sexpb, {0.1,0.2,0.3,0.4,0.5} );
    } // THEN
  } // GIVEN



} // TEST CASE

  //  for ( auto entry : bes ) { std::cout << entry << std::endl; }
