#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "discre/discre.h"
#include <iostream>

auto equal = [](auto x, auto y, double tol = 1e-5){return x == Approx(y).epsilon(tol);};

TEST_CASE( "Discrete oscillator treatment" ){
  //std::cout << std::setprecision(18);
  double temp  = 200.0012;
  GIVEN( "two oscillators" ){
    WHEN( "alpha and beta values are slightly small" ){
      std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
        alpha{0.1, 0.2, 0.4, 0.8, 1.6}, beta{0.10, 0.15, 0.30, 0.60, 1.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      double lambda_s = 2.2941534E-3, sc = 1, scaling = 1, tbeta = 2, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.9575582, 1.914953, 2.872356, 
        3.829659, 4.807700, 5.501617, 6.418366, 7.335187, 8.251386, 9.253313, 
        9.256188, 10.09717, 10.93859, 11.77685, 12.85839, 11.36789, 12.07716, 
        12.78842, 13.48618, 14.74575, 10.74546, 11.25461, 11.77105, 12.24116, 
        13.75305};

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN

    WHEN( "alpha and beta values are slightly larger" ){
      std::vector<double> osc_energies{0.035, 0.05}, osc_weights{0.2, 0.8},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }


      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tbeta = 2.0, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.7125247, 1.321226, 1.681776, 2.414681,
        3.053423, 2.268664, 3.482821, 3.839385, 5.171756, 5.720471, 2.353392, 
        3.703027, 4.256384, 6.170473, 6.777184, 1.859812, 2.977392, 3.629533, 
        5.625912, 6.256402, 1.551690, 2.499096, 3.173748, 5.098916, 5.735165};

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN

    WHEN( "oscillator weights add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.5}, osc_weights{0.4, 0.6},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tbeta = 2.0, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.8323044, 1.665653, 2.505608, 
        3.415392, 4.259058, 4.128312, 4.839682, 5.601205, 7.150784, 8.046933, 
        6.217680, 6.833610, 7.560954, 10.11139, 11.17169, 7.300605, 7.836529, 
        8.549848, 12.21146, 13.47660, 8.160909, 8.655053, 9.386204, 14.08904, 
        15.57309};

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN
    WHEN( "oscillator weights don't add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.5}, osc_weights{0.4, 0.5},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tbeta = 2.0, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);
      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.83835105, 1.6777545, 2.5238118, 
        3.4402043, 4.2899999, 4.1885133, 4.9102577, 5.6828846, 7.2550608, 
        8.1642777, 6.3563702, 6.9860384, 7.7296068, 10.336937, 11.420888, 
        7.5254485, 8.0778780, 8.8131661, 12.587555, 13.891658, 8.4646043, 
        8.9771372, 9.7354967, 14.613342, 16.15261};

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "3 or more oscillators" ){
    WHEN( "3 oscillator weights add up to 1" ){
      std::vector<double> osc_energies{0.1, 0.2, 0.3}, osc_weights{0.2, 0.3, 0.5},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tbeta = 2.0, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);
      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.8313648, 1.663253, 2.498450, 
        3.368623, 4.205679, 4.132822, 4.833326, 5.558994, 6.680081, 7.473300, 
        6.230090, 6.821939, 7.469791, 9.034025, 9.849238, 7.320823, 7.818479, 
        8.405469, 10.47719, 11.34208, 8.187767, 8.630877, 9.193504, 11.75896, 
        12.70126}; 

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN
    
    WHEN( "5 oscillators where weights don't add up to 1" ){
      std::vector<double> osc_energies{1,2,3,4,5}, osc_weights{0.5,0.4,0.3,0.2,0.1},
        alpha{2.1, 4.2, 6.4, 8.8, 10.6}, beta{1.10, 2.15, 3.30, 4.60, 5.20};
      double t_eff = 81178.935219;

      std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      double lambda_s = 2.2941534E-3, sc = 1.0, scaling = 1.0, 
             tbeta = 2.0, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);
      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
              temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.9690026, 1.938005, 2.907007, 
        3.876010, 4.845013, 5.633796, 6.572763, 7.511729, 8.450695, 9.389661, 
        9.993471, 10.90196, 11.81046, 12.71896, 13.62746, 14.02216, 14.89855, 
        15.77493, 16.65132, 17.52770, 17.91401, 18.76706, 19.62011, 20.47316, 
        21.32620}; 

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "dummy test" ){
    double sc = 1.0, scaling = 1.0,
    lambda_s = 0.26460498561058793, temp = 296.0,
    twt = 0.0, tbeta = 0.5, effectiveTemp = 572.61028482016525;
    std::vector<double> 
      alpha = {0.1, 1.0, 10},
      beta  = {0.0, 5.0, 50.0},
      oscEnergies { 0.02, 0.05 },
      oscWeights  { 0.2, 0.3 },
      sab { 6.280552E-03, 1.484229E-03, 1.781044E-29, 5.104245E-02,
            1.522726E-02, 4.744955E-20, 7.013400E-02, 8.693012E-02,
            8.708265E-10 },
      sabCorrect { 1.434008E-02, 1.499370E-03, 0.000000E+00, 6.461972E-02,
                   1.711452E-02, 1.054588E-19, 9.475565E-03, 4.652109E-02,
                   2.206891E-08 };

    auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);

    discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, oscEnergiesWeights, effectiveTemp, sab );
    THEN( "scattering law matrix is correctly changed" ){
      REQUIRE(ranges::equal(sab,sabCorrect,equal));
    } // THEN
  } // GIVEN



  GIVEN( "simple water example" ){
    double sc = 1.0, scaling = 1.0,
    lambda_s = 0.23520419644942447, temp = 296.0,
    twt = 0.055556, tbeta = 0.444444, effectiveTemp = 541.87556285322705;
    std::vector<double> 
      alpha = {0.01, 0.04, 0.08, 0.10, 0.40, 0.80, 1.00, 4.00, 8.00, 10.00},
      beta  = {0.00, 1.00, 2.00, 3.00, 4.00},
      oscEnergies { 0.205, 0.48 },
      oscWeights  { 0.166667, 0.333333 },
      sab { 1.193900E+01, 3.646870E-4, 4.920232E-4, 4.461151E-4, 1.144283E-4, 
      5.926937E+0, 1.466600E-3, 1.947206E-3, 1.763130E-3, 4.720412E-4, 
      4.152104E+00, 2.938930E-3, 3.850049E-3, 3.485422E-3, 9.678104E-4, 
      3.696742E+00, 3.674192E-3, 4.787192E-3, 4.334504E-3, 1.221669E-3, 
      1.730869E+00, 1.448399E-2, 1.787203E-2, 1.631358E-2, 5.320333E-3, 
      1.127490E+00, 3.426555E-2, 3.290205E-2, 3.048440E-2, 1.125278E-2, 
      9.694371E-1, 5.096100E-2, 3.952124E-2, 3.689802E-2, 1.431383E-2, 
      2.870226E-1, 2.015616E-1, 9.850062E-2, 9.105633E-2, 5.377092E-2, 
      1.126476E-1, 1.323797E-1, 1.095493E-1, 9.780977E-2, 7.193724E-2, 
      7.712071E-2, 1.000004E-1, 9.660049E-2, 8.905284E-2, 7.049105E-2 },
      sabCorrect {1.193441E1, 3.645467E-4, 4.918339E-4, 4.459435E-4, 
      1.143843E-4, 5.917825, 1.464345E-3, 1.944212E-3, 1.760420E-3, 4.713155E-4, 
      4.139347, 2.929900E-3, 3.838220E-3, 3.474713E-3, 9.648392E-4, 3.682550, 
      3.660087E-3, 4.768813E-3, 4.317864E-3, 1.217011E-3, 1.704442, 1.426286E-2, 
      1.759921E-2, 1.606475E-2, 5.240240E-3, 1.093324, 3.322725E-2, 3.190517E-2, 
      2.956113E-2, 1.091342E-2, 9.328565E-1, 4.903809E-2, 3.803011E-2, 
      3.550622E-2, 1.377530E-2, 2.460920E-1, 1.728181E-1, 8.445423E-2, 
      7.807188E-2, 4.610450E-2, 8.281073E-2, 9.731646E-2, 8.053320E-2, 
      7.190342E-2, 5.288448E-2, 5.249625E-2, 6.807050E-2, 6.575628E-2, 
      6.061882E-2, 4.798440E-2};

    auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);
    discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, oscEnergiesWeights, effectiveTemp, sab );
    THEN( "scattering law matrix is correctly changed" ){
      REQUIRE(ranges::equal(sab,sabCorrect,equal));
    } // THEN
  } // GIVEN







} // TEST CASE
