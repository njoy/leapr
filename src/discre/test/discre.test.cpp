#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "discre/discre.h"
#include <unsupported/Eigen/CXX11/Tensor>


void equal_vec_mega_vec( Eigen::Tensor<double,3>& sab, 
  std::vector<double> correctSab ){

  REQUIRE( sab.dimension(0)*sab.dimension(1)*sab.dimension(2) == correctSab.size() );

  int l = 0;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      for ( int k = 0; k < sab.dimension(2); ++k ){
        REQUIRE( sab(i,j,k) == Approx(correctSab[l]).epsilon(1e-5) );
	l += 1;
      }
    }
  }
}

auto populateSymSab( const std::vector<double>& alpha, const std::vector<double>& beta ){
  Eigen::Tensor<double,3> sab(alpha.size(),beta.size(),1);
  int k = 1;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      sab(i,j,0) = k;
      k += 1;
    }
  }
  return sab;
}



TEST_CASE( "discre" ){
  GIVEN( "two oscillators" ){
    WHEN( "alpha and beta values are slightly small" ){
      std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
        alpha{0.1, 0.2, 0.4, 0.8, 1.6}, beta{0.10, 0.15, 0.30, 0.60, 1.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, tev = 1.723477E-2,
             tbeta = 2.0, twt = 0.3;

      int itemp = 0;

      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.9575582, 1.914953, 2.872356, 
        3.829659, 4.807700, 5.501617, 6.418366, 7.335187, 8.251386, 9.253313, 
        9.256188, 10.09717, 10.93859, 11.77685, 12.85839, 11.36789, 12.07716, 
        12.78842, 13.48618, 14.74575, 10.74546, 11.25461, 11.77105, 12.24116, 
        13.75305};

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN

    WHEN( "alpha and beta values are slightly larger" ){
      std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;

      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.7125247, 1.321226, 1.681776, 2.414681,
        3.053423, 2.268664, 3.482821, 3.839385, 5.171756, 5.720471, 2.353392, 
        3.703027, 4.256384, 6.170473, 6.777184, 1.859812, 2.977392, 3.629533, 
        5.625912, 6.256402, 1.551690, 2.499096, 3.173748, 5.098916, 5.735165};

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN

    WHEN( "oscillator weights add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.5}, osc_weights{0.4, 0.6},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;

      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.8323044, 1.665653, 2.505608, 
        3.415392, 4.259058, 4.128312, 4.839682, 5.601205, 7.150784, 8.046933, 
        6.217680, 6.833610, 7.560954, 10.11139, 11.17169, 7.300605, 7.836529, 
        8.549848, 12.21146, 13.47660, 8.160909, 8.655053, 9.386204, 14.08904, 
        15.57309};

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN
    WHEN( "oscillator weights don't add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.5}, osc_weights{0.4, 0.5},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;

      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.83835105, 1.6777545, 2.5238118, 
        3.4402043, 4.2899999, 4.1885133, 4.9102577, 5.6828846, 7.2550608, 
        8.1642777, 6.3563702, 6.9860384, 7.7296068, 10.336937, 11.420888, 
        7.5254485, 8.0778780, 8.8131661, 12.587555, 13.891658, 8.4646043, 
        8.9771372, 9.7354967, 14.613342, 16.15261};

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "3 or more oscillators" ){
    WHEN( "3 oscillator weights add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.2, 0.3}, osc_weights{0.2, 0.3, 0.5},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;
      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.8313648, 1.663253, 2.498450, 
        3.368623, 4.205679, 4.132822, 4.833326, 5.558994, 6.680081, 7.473300, 
        6.230090, 6.821939, 7.469791, 9.034025, 9.849238, 7.320823, 7.818479, 
        8.405469, 10.47719, 11.34208, 8.187767, 8.630877, 9.193504, 11.75896, 
        12.70126}; 

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN
    
    WHEN( "5 oscillators where weights don't add up to 1" ){
      std::vector<double> osc_energies{1,2,3,4,5}, osc_weights{0.5,0.4,0.3,0.2,0.1},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20}, 
        t_eff_vec{81178.935219}, temp_vec{200.0};

      auto sym_sab = populateSymSab( alpha, beta );

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tev = 1.723477E-2, tbeta = 2.0, twt = 0.3;
      int itemp = 0;

      discre( itemp, sc, scaling, tev, lambda_s, twt, tbeta, alpha, beta, 
              temp_vec, osc_energies, osc_weights, t_eff_vec, sym_sab );

      std::vector<double> correctSymSab {0.9690026, 1.938005, 2.907007, 
        3.876010, 4.845013, 5.633796, 6.572763, 7.511729, 8.450695, 9.389661, 
        9.993471, 10.90196, 11.81046, 12.71896, 13.62746, 14.02216, 14.89855, 
        15.77493, 16.65132, 17.52770, 17.91401, 18.76706, 19.62011, 20.47316, 
        21.32620}; 

      THEN( "scattering law matrix is correctly changed" ){
        equal_vec_mega_vec( sym_sab, correctSymSab );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
