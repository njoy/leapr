#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../catch.hpp"
#include "coher.h"

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



TEST_CASE( "coher" ){

  int iel = 1, npr = 1, maxb = 60000;
  std::vector<double> b ( 60000, 0.0 );
  double emax = 5.0;

  GIVEN( "beryllium oxide is the requested material" ){
    iel = 3;
    WHEN( "1 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
	coher( iel, npr, maxb, b, emax );
	//std::vector<double> bVals { 
        //for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
	//std::cout << b[0] << "    " << b[1] << "     " << b[2] << std::endl;
      } // THEN
    } // WHEN
    WHEN( "3 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
        npr = 3;
	coher( iel, npr, maxb, b, emax );
        //std::vector<double> bVals { 
       // for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "aluminum is the requested material" ){
    iel = 4;
    WHEN( "1 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 3.759017E-3, 5.508125E-3, 5.012023E-3, 
          3.577632E-3, 1.002404E-2, 5.059536E-3, 1.378306E-2, 8.629574E-3, 
          1.503607E-2, 2.754062E-3, 2.004809E-2, 1.788816E-3, 2.380711E-2, 
          6.566121E-3, 2.506011E-2, 6.399863E-3, 3.007214E-2, 5.842249E-3, 
          3.383116E-2, 7.344167E-3, 4.009619E-2, 2.529768E-3, 4.385520E-2, 
          9.675684E-3, 4.510821E-2, 5.962720E-3, 5.012023E-2, 4.525386E-3, 
          5.387925E-2, 4.364670E-3, 5.513226E-2, 4.314787E-3, 6.014428E-2, 
          1.377031E-3, 6.390330E-2, 8.015499E-3, 6.515631E-2, 3.969026E-3, 
          7.016833E-2, 7.649299E-3 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
    WHEN( "3 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
        npr = 3;
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 3.759017E-3, 1.836041E-3, 5.012023E-3, 
          1.192544E-3, 1.002404E-2, 1.686512E-3, 1.378306E-2, 2.876524E-3, 
          1.503607E-2, 9.180209E-4, 2.004809E-2, 5.962720E-4, 2.380711E-2, 
          2.188707E-3, 2.506011E-2, 2.133287E-3, 3.007214E-2, 1.947416E-3, 
          3.383116E-2, 2.448055E-3, 4.009619E-2, 8.432560E-4, 4.385520E-2, 
          3.225228E-3, 4.510821E-2, 1.987573E-3, 5.012023E-2, 1.508462E-3, 
          5.387925E-2, 1.454890E-3, 5.513226E-2, 1.438262E-3, 6.014428E-2, 
          4.590104E-4, 6.390330E-2, 2.671833E-3, 6.515631E-2, 1.323008E-3, 
          7.016833E-2, 2.549766E-3 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
  } // GIVEN



  GIVEN( "lead is the requested material" ){
    iel = 5;
    WHEN( "1 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 2.514103E-3, 2.464174E-3, 3.352138E-3, 
          1.600528E-3, 6.704277E-3, 2.263488E-3, 9.218381E-3, 3.860619E-3, 
          1.005641E-2, 1.232087E-3, 1.340855E-2, 8.002641E-4, 1.592265E-2, 
          2.937491E-3, 1.676069E-2, 2.863112E-3, 2.011283E-2, 2.613651E-3, 
          2.262693E-2, 3.285566E-3, 2.681710E-2, 1.131744E-3, 2.933121E-2, 
          4.328618E-3, 3.016924E-2, 2.667547E-3, 3.352138E-2, 2.024525E-3, 
          3.603549E-2, 1.952626E-3, 3.687352E-2, 1.930309E-3, 4.022566E-2, 
          6.160436E-4, 4.273976E-2, 3.585900E-3, 4.357780E-2, 1.775626E-3, 
          4.692994E-2, 3.422073E-3 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
    WHEN( "5 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
        npr = 5;
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 2.514103E-3, 4.928349E-4, 3.352138E-3, 
          3.201056E-4, 6.704277E-3, 4.526977E-4, 9.218381E-3, 7.721239E-4, 
          1.005641E-2, 2.464174E-4, 1.340855E-2, 1.600528E-4, 1.592265E-2, 
          5.874982E-4, 1.676069E-2, 5.726224E-4, 2.011283E-2, 5.227303E-4, 
          2.262693E-2, 6.571132E-4, 2.681710E-2, 2.263488E-4, 2.933121E-2, 
          8.657237E-4, 3.016924E-2, 5.335094E-4, 3.352138E-2, 4.049051E-4, 
          3.603549E-2, 3.905252E-4, 3.687352E-2, 3.860619E-4, 4.022566E-2, 
          1.232087E-4, 4.273976E-2, 7.171801E-4, 4.357780E-2, 3.551253E-4, 
          4.692994E-2, 6.844146E-4 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "iron is the requested material" ){
    iel = 6;
    WHEN( "1 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 5.000504E-3, 8.711438E-2, 1.000100E-2, 
          3.079958E-2, 1.500151E-2, 0.100591032, 2.000201E-2, 4.355719E-2, 
          2.500252E-2, 7.791747E-2, 3.000302E-2, 2.370953E-2, 3.500353E-2, 
          0.131704577, 4.000403E-2, 1.539979E-2, 4.500453E-2, 8.711438E-2, 
          5.000504E-2, 5.509597E-2, 5.500554E-2, 5.253195E-2, 6.000605E-2, 
          5.029551E-2, 6.500655E-2, 0.144967107, 7.500756E-2, 8.997135E-2, 
          8.000806E-2, 2.177859E-2, 8.500857E-2, 8.451337E-2, 9.000907E-2, 
          5.133264E-2, 9.500958E-2, 0.119912469, 0.100010086, 3.895873E-2, 
          0.105010591, 7.603967E-2 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
    WHEN( "5 principal scattering atoms in compound" ){
      THEN( "bragg edges vector is correctly output" ){
        npr = 5;
	coher( iel, npr, maxb, b, emax );
        std::vector<double> bVals { 5.000504E-3, 1.742287E-2, 1.000100E-2, 
          6.159917E-3, 1.500151E-2, 2.011820E-2, 2.000201E-2, 8.711438E-3, 
          2.500252E-2, 1.558349E-2, 3.000302E-2, 4.741906E-3, 3.500353E-2, 
          2.634091E-2, 4.000403E-2, 3.079958E-3, 4.500453E-2, 1.742287E-2, 
          5.000504E-2, 1.101919E-2, 5.500554E-2, 1.050639E-2, 6.000605E-2, 
          1.005910E-2, 6.500655E-2, 2.899342E-2, 7.500756E-2, 1.799427E-2, 
          8.000806E-2, 4.355719E-3, 8.500857E-2, 1.690267E-2, 9.000907E-2, 
          1.026652E-2, 9.500958E-2, 2.398249E-2, 0.100010086, 7.791747E-3, 
          0.105010591, 1.520793E-2 };
        for ( auto i = 0; i < bVals.size(); ++i ){ equal( b[i], bVals[i] ); }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
