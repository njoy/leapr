#include "catch.hpp"
#include "simple/continuous/sumOverTn.h"



TEST_CASE( "sumOverTn" ){
  /*
  GIVEN( "Only one iteration is performed (only want T0)" ){
    THEN( "We only see the effects of the T0(b) = delta(b) term" ){
      {
        std::vector<double> alphas {2,4}, betas {0,1,2,3,4}, T1 {1,2,3,2,1};
        double lambda_s = 0.5;
        int N = 1;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( exp(-1) == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( exp(-2) == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        for ( size_t b = 1; b < betas.size(); ++b ){
          REQUIRE( 0.0 == Approx(sab[0*betas.size()+b]).epsilon(1e-6) );
          REQUIRE( 0.0 == Approx(sab[1*betas.size()+b]).epsilon(1e-6) );
        }
      }
      {
        std::vector<double> alphas {0.1,0.4}, betas {0,2}, T1 {1,3};
        double lambda_s = 0.75;
        int N = 1;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( 0.92774348632 == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.74081822068 == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.0 == Approx(sab[0*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 0.0 == Approx(sab[1*betas.size()+1]).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "Two iterations are performed (T0 and T1)" ){
    THEN( "T0(b) = delta(b) and T1 which is already defined are used" ){
      {
        std::vector<double> alphas {2,4}, betas {0,1,2,3,4}, T1 {1,2,3,2,1};
        double lambda_s = 0.5;
        int N = 2;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( 2.0*exp(-1) == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 2.0*exp(-1) == Approx(sab[0*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 3.0*exp(-1) == Approx(sab[0*betas.size()+2]).epsilon(1e-6) );
        REQUIRE( 2.0*exp(-1) == Approx(sab[0*betas.size()+3]).epsilon(1e-6) );
        REQUIRE( 1.0*exp(-1) == Approx(sab[0*betas.size()+4]).epsilon(1e-6) );

        REQUIRE( 3.0*exp(-2) == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 4.0*exp(-2) == Approx(sab[1*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 6.0*exp(-2) == Approx(sab[1*betas.size()+2]).epsilon(1e-6) );
        REQUIRE( 4.0*exp(-2) == Approx(sab[1*betas.size()+3]).epsilon(1e-6) );
        REQUIRE( 2.0*exp(-2) == Approx(sab[1*betas.size()+4]).epsilon(1e-6) );
      }
      {
        std::vector<double> alphas {0.1,0.4}, betas {0,2}, T1 {1,3};
        double lambda_s = 0.75;
        int N = 2;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( 0.99732424780 == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.96306368688 == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.20874228442 == Approx(sab[0*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 0.66673639861 == Approx(sab[1*betas.size()+1]).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "Three iterations are performed (T0, T1, T2)" ){
    THEN( "" ){
      {
        std::vector<double> alphas {0.1,0.4}, betas {0,2}, T1 {1,3};
        double lambda_s = 0.75;
        int N = 3;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( 1.3495847061217 == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 5.4636347140466 == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.2322257914216 == Approx(sab[0*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 0.9667677779896 == Approx(sab[1*betas.size()+1]).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "Four iterations are performed (T0, T1, T2, T3)" ){
    THEN( "" ){
      {
        std::vector<double> alphas {0.1,0.4}, betas {0,2}, T1 {1,3};
        double lambda_s = 0.75;
        int N = 4;
        auto sab = sumOverTn(alphas,betas,T1,lambda_s,N);
        REQUIRE( 1.4192540142165 == Approx(sab[0*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 9.0240873517047 == Approx(sab[1*betas.size()+0]).epsilon(1e-6) );
        REQUIRE( 0.3122671781088 == Approx(sab[0*betas.size()+1]).epsilon(1e-6) );
        REQUIRE( 5.0572858862615 == Approx(sab[1*betas.size()+1]).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  */

} // TEST



