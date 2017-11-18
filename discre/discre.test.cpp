#define CATCH_CONFIG_MAIN
#include "../catch.hpp" 
#include "discre.h"
#include <iostream> 


void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec_mega_vec( std::vector<std::vector<std::vector<double>>> a, 
  std::vector<double> b ){
  REQUIRE( a.size()*a[0].size()*a[0][0].size() == b.size() );
  int i = 0;
  for ( auto a1 : a ){
    for ( auto a2 : a1 ){
      for ( auto a3 : a2 ){
        equal( a3, b[i] );
        i += 1;
      }
    }
  }
}




TEST_CASE( "discre" ){
  GIVEN( "inputs1" ){
    std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
      alpha{0.1, 0.2, 0.4, 0.8, 1.6}, beta{0.10, 0.15, 0.30, 0.60, 1.20}, 
      t_eff_vec{81178.935219}, temp_vec{200.0};

    std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
      std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));
    int i = 1;
    for ( auto a = 0; a < alpha.size(); ++a ){
      for ( auto b = 0; b < beta.size(); ++b ){
        sym_sab[a][b][0] = i;
        i += 1;
      }
    }
    double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, tev = 1.723477E-2,
           tbeta = 2.0, twt = 0.3;
    int itemp = 0;

    discre( sc, scaling, alpha, beta, tev, lambda_s, osc_energies, osc_weights, 
        tbeta, t_eff_vec, temp_vec, itemp, sym_sab, twt );

    std::vector<double> correctSymSab {0.9575582, 1.914953, 2.872356, 3.829659,
      4.807700, 5.501617, 6.418366, 7.335187, 8.251386, 9.253313, 9.256188, 
      10.09717, 10.93859, 11.77685, 12.85839, 11.36789, 12.07716, 12.78842, 
      13.48618, 14.74575, 10.74546, 11.25461, 11.77105, 12.24116, 13.75305};

    THEN( "scattering law matrix is correctly changed" ){
      equal_vec_mega_vec( sym_sab, correctSymSab );
    } // THEN
  } // GIVEN

  GIVEN( "large alpha and beta values" ){
    std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
      alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
      t_eff_vec{81178.935219}, temp_vec{200.0};

    std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
      std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));
    int i = 1;
    for ( auto a = 0; a < alpha.size(); ++a ){
      for ( auto b = 0; b < beta.size(); ++b ){
        sym_sab[a][b][0] = i;
        i += 1;
      }
    }
    double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, tev = 1.723477E-2,
           tbeta = 2.0, twt = 0.3;
    int itemp = 0;

    discre( sc, scaling, alpha, beta, tev, lambda_s, osc_energies, osc_weights, 
        tbeta, t_eff_vec, temp_vec, itemp, sym_sab, twt );

    std::vector<double> correctSymSab {0.7125247, 1.321226, 1.681776, 2.414681,
      3.053423, 2.268664, 3.482821, 3.839385, 5.171756, 5.720471, 2.353392, 
      3.703027, 4.256384, 6.170473, 6.777184, 1.859812, 2.977392, 3.629533, 
      5.625912, 6.256402, 1.551690, 2.499096, 3.173748, 5.098916, 5.735165};

    THEN( "scattering law matrix is correctly changed" ){
      equal_vec_mega_vec( sym_sab, correctSymSab );
    } // THEN
  } // GIVEN

  GIVEN( "more oscillators" ){
    WHEN( "oscillator weights don't add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.5, 1}, osc_weights{0.2, 0.7, 0.1},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
      std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));

      int i = 1;
      for ( auto a = 0; a < alpha.size(); ++a ){
        for ( auto b = 0; b < beta.size(); ++b ){
          sym_sab[a][b][0] = i;
          i += 1;
        }
      }
      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;

      discre( sc, scaling, alpha, beta, tev, lambda_s, osc_energies, 
          osc_weights, tbeta, t_eff_vec, temp_vec, itemp, sym_sab, twt );

      std::vector<double> correctSymSab {0.7125247, 1.321226, 1.681776, 2.414681, 3.053423, 2.268664, 3.482821, 3.839385, 5.171756, 5.720471, 2.353392, 3.703027, 4.256384, 6.170473, 6.777184, 1.859812, 2.977392, 3.629533, 5.625912, 6.256402, 1.551690, 2.499096, 3.173748, 5.098916, 5.735165};

      //THEN( "scattering law matrix is correctly changed" ){
      //  equal_vec_mega_vec( sym_sab, correctSymSab );
     // } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
