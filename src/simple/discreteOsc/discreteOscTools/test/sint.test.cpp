#include "catch.hpp"
#include "simple/discreteOsc/discreteOscTools/sint.h"

TEST_CASE( "sint" ){
  GIVEN( "inputs" ){
    double x, alpha, wt = 2.3, tbart = 407.4545311;
    int b = 4, nbx = 10;
    std::vector<double> bex {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0};
    std::vector<double> rdbex  {1.66666667, 3.33333333, 6.666666667, 20.0, 5.0,
      20.0, 6.66666667, 3.33333333, 1.66666667, 0.0, 0.0};
    std::vector<double> betan {0.1, 0.15, 0.3, 0.6, 1.2};
    std::vector<double> sex  {5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.7214159, 
      2.2224546, 2.1952465, 1.5059710, 0.0};
  
    double sintOut;


    WHEN( "beta value (larger than grid range), and negative alpha" ){
      THEN( "0.0 is returned because -alpha is invalid for SCT approximation" ){
        x = -1.201; alpha = -0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(0.0).epsilon(1e-6) );
        x = 1.201; alpha = -1e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "positive beta (larger than grid range), and positive alpha" ){
      THEN( "SCT approx. made, treating beta as -beta, then mult by e^-beta" ){
        x = 1.201; alpha = 0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(2.548608E-4).epsilon(1e-6) );
      } // THEN
    } // WHEN


    WHEN( "negtive beta (larger than grid range), and positive alpha " ){
      THEN( "SCT approx. made" ){
        x = -1.201; alpha = 1e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(1.654202E-16).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "Desired beta value is within range and is a bisection point" ){
      THEN( "Corresponding tabulated value is returned" ){
        x = -0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(1.0).epsilon(1e-6) );
        x = -0.15;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(2.0).epsilon(1e-6) );
        x = -0.3;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(3.0).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "Different alpha values provided for same interpolation problem" ){
      THEN( "no change in output since alpha only important in SCT approx." ){ 
        x = -0.55; alpha = 0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(3.812737215).epsilon(1e-6) );
        alpha = 1.0e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(3.812737215).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "bisect left to find desired point between two positive values" ){
      THEN( "interpolate between the two to get correct value" ){
        x = -0.35;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(3.14734517).epsilon(1e-6) );

        x = -0.65; 
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(4.07507702).epsilon(1e-6) );

      } // THEN
    } // WHEN
     
    WHEN( "bisect right to find desired point between two positive values" ){
      THEN( "interpolate between the two to get correct value" ){
        x = -0.11; alpha = 1.0e-5; 
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(1.148698345).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "bisect right to find desired point between a neg val and pos val" ){
      THEN( "log(negative value to the left of desired point) is set to 0.0" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};
        x = -0.11; alpha = 1.0e-5; 
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(1.921947E-98).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "bisect and interpolate desired point, but resulting exp. is small" ){
      THEN( "output value of zero" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};

        rdbex = {1.66666667, 3.33333333, 6.666666667, 25.0, 5.0, 25.0, 
          6.66666667, 3.33333333, 1.66666667, 0.0, 0.0};
        x = -0.11; alpha = 1.0e-5; 
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        REQUIRE( sintOut == Approx(0.0).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
