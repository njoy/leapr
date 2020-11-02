#include "catch.hpp"
#include "translational/translational_util/s_table_generation.h"
#include "generalTools/testing.h"
#include <range/v3/all.hpp>

TEST_CASE( "diffusion s-table generation" ){
  GIVEN( "trans. weight, alpha, diffusion const, delta and empty vectors" ){
    THEN( "the beta and translational scattering vectors are filled" ){
      double c = 2.0, delta = 0.2, trans_weight = 1.5, alpha = 0.1;
      int ndmax = 100, nsd;
      std::vector<double> sd(ndmax,0.0);
      nsd = getDiffusion( trans_weight, alpha, ndmax, delta, sd, c );
      std::vector<double> correct_sd {0.8939586849, 0.8460637169, 0.6261350235, 
      0.4069039776, 0.2516455797, 0.1539260520, 9.45900071E-2, 5.87182586E-2, 
      3.68702007E-2, 2.34112959E-2, 1.50191406E-2, 9.72506689E-3, 6.34946619E-3, 
      4.17624927E-3, 2.76498632E-3, 1.84142429E-3, 1.23283758E-3, 8.29317113E-4, 
      5.60273054E-4, 3.79988371E-4, 2.58631358E-4, 1.76603655E-4, 1.20950854E-4, 
      8.30626003E-5, 5.71870387E-5, 3.94642033E-5, 2.72928875E-5, 1.89134301E-5, 
      1.31312649E-5, 9.13281730E-6, 6.36233337E-6, 4.43912567E-6, 3.10175614E-6, 
      2.17025290E-6, 1.52044914E-6, 1.06650140E-6, 7.48946055E-7, 5.26517742E-7, 
      3.70531078E-7, 2.61013523E-7, 1.84037923E-7, 1.29878799E-7, 9.17355875E-8, 
      6.48469122E-8, 4.58750947E-8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      REQUIRE( nsd == 45 );
      REQUIRE(ranges::equal(sd,correct_sd,equal));
    } // THEN
  } // GIVEN
  GIVEN( "translational weight, scaled alpha, spacing and empty vectors" ){
    double c = 0.3, delta = 1.8, trans_weight = 0.5, alpha = 0.8;
    int ndmax = 20, nsd;
    std::vector<double> sd(ndmax,0.0);
    THEN( "the beta and translational scattering vectors are filled" ){
      std::vector<double> correct_sd{1.3891497, 3.5359670E-2, 9.8039660E-3, 
      4.4150665E-3, 2.4152380E-3, 1.4673486E-3, 9.5188023E-4, 6.4581826E-4, 
      4.5268825E-4, 3.2527107E-4, 2.3830869E-4, 1.7735539E-4, 1.3370720E-4, 
      1.0189672E-4, 7.8370923E-5, 6.0754890E-5, 4.7423166E-5, 3.7240618E-5, 
      2.9400700E-5, 0};
      nsd = getDiffusion( trans_weight, alpha, ndmax, delta, sd, c );
      REQUIRE( nsd == 19 );
      REQUIRE(ranges::equal(sd,correct_sd,equal));
    } // THEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "free gas s-table generation" ){
  GIVEN( "trans. weight, alpha, delta and empty vectors" ){
    THEN( "the beta and translational scattering vectors are filled" ){
      double delta = 1.4, trans_weight = 0.5, alpha = 0.8;
      int ndmax = 10, nsd;
      std::vector<double> sd(ndmax,0.0);
      nsd = getFreeGas( trans_weight, alpha, ndmax, delta, sd );
      std::vector<double> correct_sd {0.40358556, 0.23874321, 1.2187230E-2, 
        5.3685572E-5, 2.0407449E-8, 0.0, 0.0, 0.0, 0.0, 0.0};
      REQUIRE( nsd == 5 );
      REQUIRE(ranges::equal(sd,correct_sd,equal));
    } // THEN
  } // GIVEN
  GIVEN( "translational weight, scaled alpha, spacing and empty vectors" ){
    double delta = 1.2, trans_weight = 0.3, alpha = 1.2;
    int ndmax = 20, nsd;
    std::vector<double> sd(ndmax,0.0);
    THEN( "the beta and translational scattering vectors are filled" ){
      std::vector<double> correct_sd{0.429692, 0.2880312, 2.6129596687E-2, 
        3.20802E-4, 5.33031E-7, 1.1986111752398333E-10, 3.647675E-15, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      
      nsd = getFreeGas( trans_weight, alpha, ndmax, delta, sd );
      REQUIRE( nsd == 7 );
      REQUIRE(ranges::equal(sd,correct_sd,equal));
    } // THEN
  } // GIVEN
} // TEST CASE

           
