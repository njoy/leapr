#include "catch.hpp"
#include "trans/trans_util/s_table_generation.h"

void equal1( double a, double b ){
    if( b == 0 ){ REQUIRE( (a-b) < 1e-6 ); }
    if( b != 0 ){ REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 ); }
}

void equal1_vec( std::vector<double> a, std::vector<double> b ){
    REQUIRE( a.size() == b.size() );
    for ( size_t i = 0; i < a.size(); ++i ){
        equal1( a[i], b[i] );
    }
}



TEST_CASE( "diffusion s-table generation" ){
  GIVEN( "trans. weight, alpha, diffusion const, delta and empty vectors" ){
    THEN( "the beta and translational scattering vectors are filled" ){
      double c = 2.0, delta = 0.2, trans_weight = 1.5, alpha = 0.1;
      int ndmax = 100, nsd;
      std::vector<double> sd(ndmax,0.0);
      nsd = diffusion_s_table( trans_weight, alpha, ndmax, delta, sd, c );
      std::vector<double> correct_sd {0.8939586849, 0.8460637169, 
        0.6261350235,    0.4069039776,    0.2516455797,    0.1539260520, 
        9.4590007147E-2, 5.8718258693E-2, 3.6870200714E-2, 2.3411295948E-2, 
        1.5019140685E-2, 9.7250668921E-3, 6.3494661925E-3, 4.1762492756E-3, 
        2.7649863261E-3, 1.8414242905E-3, 1.2328375818E-3, 8.2931711329E-4, 
        5.6027305464E-4, 3.7998837198E-4, 2.5863135811E-4, 1.7660365598E-4, 
        1.2095085418E-4, 8.3062600333E-5, 5.7187038779E-5, 3.9464203353E-5, 
        2.7292887506E-5, 1.8913430113E-5, 1.3131264933E-5, 9.1328173051E-6, 
        6.3623333797E-6, 4.4391256752E-6, 3.1017561487E-6, 2.1702529083E-6, 
        1.5204491464E-6, 1.0665014088E-6, 7.4894605566E-7, 5.2651774286E-7, 
        3.7053107838E-7, 2.6101352369E-7, 1.8403792387E-7, 1.2987879955E-7, 
        9.1735587540E-8, 6.4846912248E-8, 4.5875094732E-8, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      REQUIRE( nsd == 45 );
      equal1_vec(sd,correct_sd);
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
        2.9400700E-5, 0.0};
      nsd = diffusion_s_table( trans_weight, alpha, ndmax, delta, sd, c );
      equal1( nsd, 19 );
      equal1_vec(sd,correct_sd);

    } // THEN
  } // GIVEN
} // TEST CASE

 

TEST_CASE( "free gas s-table generation" ){
  GIVEN( "trans. weight, alpha, delta and empty vectors" ){
    THEN( "the beta and translational scattering vectors are filled" ){
      double delta = 1.4, trans_weight = 0.5, alpha = 0.8;
      int ndmax = 10, nsd;
      std::vector<double> sd(ndmax,0.0);
      nsd = free_gas_s_table( trans_weight, alpha, ndmax, delta, sd );
      std::vector<double> correct_sd {0.40358556, 0.23874321, 1.2187230E-2, 
        5.3685572E-5, 2.0407449E-8, 0.0, 0.0, 0.0, 0.0, 0.0};
      equal1( nsd, 5 );
      equal1_vec(sd,correct_sd);
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
      
      nsd = free_gas_s_table( trans_weight, alpha, ndmax, delta, sd );
      equal1( nsd, 7 );
      equal1_vec(sd,correct_sd);

    } // THEN
  } // GIVEN
} // TEST CASE

           
