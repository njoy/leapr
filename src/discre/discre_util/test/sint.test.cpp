#include "catch.hpp"
#include "discre/discre_util/sint.h"

TEST_CASE( "sint" ){
  GIVEN( "inputs" ){
    double x, alpha, wt = 2.3, tbart = 407.4545311, sintOut;
    int b = 4, nbx = 10;
    std::vector<double> 
      bex   {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0},
      rdbex {1.666667,3.333333,6.666667,20,5,20,6.666667,3.33333,1.666667,0,0},
      betan {0.1, 0.15, 0.3, 0.6, 1.2},
      sex   {5, 4, 3, 2, 1, 1, 1.7214159, 2.2224546, 2.1952465, 1.505971, 0};

    WHEN( "beta value (larger than grid range), and negative alpha" ){
      THEN( "0.0 is returned because -alpha is invalid for SCT approximation" ){
        sintOut = sint(x=-1.201,bex,rdbex,sex,betan,b,alpha=-0.1*wt,tbart,nbx);
        REQUIRE( 0.0 == Approx(sintOut).epsilon(1e-6) );

        sintOut = sint(x=1.201,bex,rdbex,sex,betan,b,alpha=-1e-5*wt,tbart,nbx);
        REQUIRE( 0.0 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "positive beta (larger than grid range), and positive alpha" ){
      THEN( "SCT approx. made, treating beta as -beta, then mult by e^-beta" ){
        sintOut = sint(x=1.201,bex,rdbex,sex,betan,b,alpha=0.1*wt,tbart,nbx);
        REQUIRE( 2.548608E-4 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "negtive beta (larger than grid range), and positive alpha " ){
      THEN( "SCT approx. made" ){
        sintOut = sint(x=-1.201,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 1.654202E-16 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "Desired beta value is within range and is a bisection point" ){
      THEN( "Corresponding tabulated value is returned" ){
        sintOut = sint(x=-0.1,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 1.0 == Approx(sintOut).epsilon(1e-6) );

        sintOut = sint(x=0-.15,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 2.0 == Approx(sintOut).epsilon(1e-6) );

        sintOut = sint(x=-0.3,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 3.0 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN
    
    WHEN( "Different alpha values provided for same interpolation problem" ){
      THEN( "no change in output since alpha only important in SCT approx." ){ 
        sintOut = sint(x=-0.55,bex,rdbex,sex,betan,b,alpha=0.1*wt,tbart,nbx);
        REQUIRE( 3.812737215 == Approx(sintOut).epsilon(1e-6) );
        sintOut = sint(x=-0.55,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 3.812737215 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "bisect left to find desired point between two positive values" ){
      THEN( "interpolate between the two to get correct value" ){
        sintOut = sint(x=-0.35,bex,rdbex,sex,betan,b,alpha=0.1*wt,tbart,nbx);
        REQUIRE( 3.14734517 == Approx(sintOut).epsilon(1e-6) );

        sintOut = sint(x=-0.65,bex,rdbex,sex,betan,b,alpha=0.1*wt,tbart,nbx);
        REQUIRE( 4.07507702 == Approx(sintOut).epsilon(1e-6) );

      } // THEN
    } // WHEN
     
    WHEN( "bisect right to find desired point between two positive values" ){
      THEN( "interpolate between the two to get correct value" ){
        sintOut = sint(x=-0.11,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 1.148698345 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "bisect right to find desired point between a neg val and pos val" ){
      THEN( "log(negative value to the left of desired point) is set to 0.0" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};
        sintOut = sint(x=-0.11,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 1.921947E-98 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "bisect and interpolate desired point, but resulting exp. is small" ){
      THEN( "output value of zero" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};

        rdbex = {1.666667,3.33333,6.666667,25,5,25,6.666667,3.33333,1.666667,0,0};
        sintOut = sint(x=-0.11,bex,rdbex,sex,betan,b,alpha=1e-5*wt,tbart,nbx);
        REQUIRE( 0.0 == Approx(sintOut).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
