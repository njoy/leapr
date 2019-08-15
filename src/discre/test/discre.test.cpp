#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "discre/discre.h"
#include <iostream>

auto equal = [](auto x, auto y, double tol = 1e-5){return x == Approx(y).epsilon(tol);};

TEST_CASE( "Discrete oscillator treatment" ){
  GIVEN( "Test material" ){
    double sc = 1.0, scaling = 1.0,
    lambda_s = 0.26460498561058793, temp = 296.0,
    twt = 0.0, tbeta = 0.5, effectiveTemp = 572.61028482016525;

    WHEN( "Single oscillator is present" ){
      std::vector<double> oscEnergies {0.3}, oscWeights {0.5};
      auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);
      AND_WHEN( "Few alpha, beta values spanning 0.1 -> 10's" ){
        //std::vector<double> 
          //alpha = {0.1, 1.0, 10},
          //beta  = {0.0, 5.0, 50.0},
          //sab { 6.28055E-3, 1.48422E-3, 1.78104E-29, 5.10424E-2, 1.52272E-2, 
          //      4.74495E-20, 7.01340E-2, 8.69301E-2, 8.7082E-10 },
          //correct  { 1.43400E-2, 1.49937E-3, 0.00000E+00, 6.46197E-2,
          //           1.71145E-2, 1.05458E-19, 9.47556E-3, 4.65210E-2,
          //           2.20689E-8 };
          //discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
          //        oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
           //REQUIRE(ranges::equal(sab,correct,equal));
        } // THEN
      } // AND WHEN
      AND_WHEN( "Very small alpha, beta values" ){
       // std::vector<double> 
       // alpha = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2},
       // beta  = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2},
       // sab { 6.42751E-7, 6.42911E-7, 6.42780E-7, 6.42911E-7, 6.43074E-7, 
       // 6.44381E-7, 6.46015E-7, 6.59084E-7, 3.21339E-5, 3.21419E-5, 3.21354E-5, 
       // 3.21419E-5, 3.21501E-5, 3.22154E-5, 3.22971E-5, 3.29504E-5, 6.42738E-6, 
       // 6.42898E-6, 6.42767E-6, 6.42898E-6, 6.43061E-6, 6.44368E-6, 6.46001E-6, 
       // 6.59070E-6, 3.21339E-5, 3.21419E-5, 3.21354E-5, 3.21419E-5, 3.21501E-5, 
       // 3.22154E-5, 3.22971E-5, 3.29504E-5, 6.42604E-5, 6.42764E-5, 6.42633E-5, 
       // 6.42764E-5, 6.42927E-5, 6.44234E-5, 6.45867E-5, 6.58933E-5, 3.21004E-4, 
       // 3.21084E-4, 3.21019E-4, 3.21084E-4, 3.21166E-4, 3.21819E-4, 3.22634E-4, 
       // 3.29161E-4, 6.41267E-4, 6.41427E-4, 6.41296E-4, 6.41427E-4, 6.41590E-4, 
       // 6.42893E-4, 6.44523E-4, 6.57560E-4, 3.17680E-3, 3.17760E-3, 3.17695E-3, 
       // 3.17760E-3, 3.17840E-3, 3.18486E-3, 3.19292E-3, 3.25746E-3 },
       // correct  { 6.42745E-7, 6.42905E-7, 6.42775E-7, 6.42905E-7, 6.43069E-7, 
       // 6.44376E-7, 6.46009E-7, 6.59078E-7, 3.21196E-5, 3.21276E-5, 3.21211E-5, 
       // 3.21276E-5, 3.21358E-5, 3.22011E-5, 3.22827E-5, 3.29358E-5, 6.42681E-6, 
       // 6.42841E-6, 6.42710E-6, 6.42841E-6, 6.43004E-6, 6.44311E-6, 6.45944E-6, 
       // 6.59012E-6, 3.21196E-5, 3.21276E-5, 3.21211E-5, 3.21276E-5, 3.21358E-5, 
       // 3.22011E-5, 3.22827E-5, 3.29358E-5, 6.42034E-5, 6.42194E-5, 6.42064E-5, 
       // 6.42194E-5, 6.42357E-5, 6.43663E-5, 6.45294E-5, 6.58349E-5, 3.19784E-4, 
       // 3.19863E-4, 3.19798E-4, 3.19863E-4, 3.19945E-4, 3.20597E-4, 3.21415E-4, 
       // 3.28138E-4, 7.14086E-4, 7.14265E-4, 7.14119E-4, 7.14265E-4, 7.14451E-4, 
       // 7.16049E-4, 7.18329E-4, 7.48380E-4, 7.16740E-3, 7.16917E-3, 7.16773E-3, 
       // 7.16917E-3, 7.17098E-3, 7.18566E-3, 7.20446E-3, 7.37363E-3 };
       // discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
       //         oscEnergiesWeights, effectiveTemp, sab );
       // THEN( "scattering law matrix is correctly changed" ){
       //   REQUIRE(ranges::equal(sab,correct,equal));
       // } // THEN
      } // AND WHEN
    } // WHEN

    WHEN( "Two oscillators present" ){
      std::vector<double> oscEnergies {0.02, 0.05}, oscWeights {0.2, 0.3};
      auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);
      AND_WHEN( "Few alpha, beta values spanning 0.1 -> 10's" ){
        std::vector<double> 
          alpha = {0.1, 1.0, 10},
          beta  = {0.0, 5.0, 50.0},
          sab { 6.28055E-3, 1.48422E-3, 1.78104E-29, 5.10424E-2, 1.52272E-2, 
                4.74495E-20, 7.01340E-2, 8.69301E-2, 8.7082E-10 },
          correct  { 1.43400E-2, 1.49937E-3, 0.00000E+00, 6.46197E-2,
                     1.71145E-2, 1.05458E-19, 9.47556E-3, 4.65210E-2,
                     2.20689E-8 };
          discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                  oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
          REQUIRE(ranges::equal(sab,correct,equal));
          } // THEN
      } // AND WHEN
      AND_WHEN( "Very small alpha, beta values" ){
        std::vector<double> 
        alpha = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2},
        beta  = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2},
        sab { 6.42751E-7, 6.42911E-7, 6.42780E-7, 6.42911E-7, 6.43074E-7, 
        6.44381E-7, 6.46015E-7, 6.59084E-7, 3.21339E-5, 3.21419E-5, 3.21354E-5, 
        3.21419E-5, 3.21501E-5, 3.22154E-5, 3.22971E-5, 3.29504E-5, 6.42738E-6, 
        6.42898E-6, 6.42767E-6, 6.42898E-6, 6.43061E-6, 6.44368E-6, 6.46001E-6, 
        6.59070E-6, 3.21339E-5, 3.21419E-5, 3.21354E-5, 3.21419E-5, 3.21501E-5, 
        3.22154E-5, 3.22971E-5, 3.29504E-5, 6.42604E-5, 6.42764E-5, 6.42633E-5, 
        6.42764E-5, 6.42927E-5, 6.44234E-5, 6.45867E-5, 6.58933E-5, 3.21004E-4, 
        3.21084E-4, 3.21019E-4, 3.21084E-4, 3.21166E-4, 3.21819E-4, 3.22634E-4, 
        3.29161E-4, 6.41267E-4, 6.41427E-4, 6.41296E-4, 6.41427E-4, 6.41590E-4, 
        6.42893E-4, 6.44523E-4, 6.57560E-4, 3.17680E-3, 3.17760E-3, 3.17695E-3, 
        3.17760E-3, 3.17840E-3, 3.18486E-3, 3.19292E-3, 3.25746E-3 },
        correct  { 6.42745E-7, 6.42905E-7, 6.42775E-7, 6.42905E-7, 6.43069E-7, 
        6.44376E-7, 6.46009E-7, 6.59078E-7, 3.21196E-5, 3.21276E-5, 3.21211E-5, 
        3.21276E-5, 3.21358E-5, 3.22011E-5, 3.22827E-5, 3.29358E-5, 6.42681E-6, 
        6.42841E-6, 6.42710E-6, 6.42841E-6, 6.43004E-6, 6.44311E-6, 6.45944E-6, 
        6.59012E-6, 3.21196E-5, 3.21276E-5, 3.21211E-5, 3.21276E-5, 3.21358E-5, 
        3.22011E-5, 3.22827E-5, 3.29358E-5, 6.42034E-5, 6.42194E-5, 6.42064E-5, 
        6.42194E-5, 6.42357E-5, 6.43663E-5, 6.45294E-5, 6.58349E-5, 3.19784E-4, 
        3.19863E-4, 3.19798E-4, 3.19863E-4, 3.19945E-4, 3.20597E-4, 3.21415E-4, 
        3.28138E-4, 7.14086E-4, 7.14265E-4, 7.14119E-4, 7.14265E-4, 7.14451E-4, 
        7.16049E-4, 7.18329E-4, 7.48380E-4, 7.16740E-3, 7.16917E-3, 7.16773E-3, 
        7.16917E-3, 7.17098E-3, 7.18566E-3, 7.20446E-3, 7.37363E-3 };
        discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
          REQUIRE(ranges::equal(sab,correct,equal));
        } // THEN
      } // AND WHEN
      /*
      AND_WHEN( "Rather large alpha, beta values" ){
        std::vector<double> 
        alpha = {1, 10, 20, 40, 60, 80, 100, 120},
        beta  = {0, 10, 20, 40, 60, 80, 100, 120},
        sab { 5.104245E-02, 4.175489E-04, 1.214367E-07, 1.074175E-15,
1.466177E-24, 0.000000E+00, 0.000000E+00, 0.000000E+00,
7.013400E-02, 3.705110E-02, 1.691192E-03, 1.947028E-07,
2.578293E-12, 8.559867E-18, 1.028888E-23, 2.333784E-30,
1.962248E-02, 6.267878E-02, 1.627107E-02, 4.809242E-05,
1.391085E-08, 9.202806E-13, 2.072519E-17, 2.001354E-22,
1.613300E-03, 2.783514E-02, 4.488554E-02, 4.289867E-03,
3.370974E-05, 5.423767E-08, 2.750372E-11, 5.645936E-15,
1.708965E-04, 5.967059E-03, 2.728335E-02, 2.154630E-02,
1.298896E-03, 1.502281E-05, 5.211192E-08, 7.040755E-11,
1.982916E-05, 1.000116E-03, 9.042812E-03, 3.190938E-02,
8.482092E-03, 4.165801E-04, 5.932873E-06, 3.200259E-08,
2.392988E-06, 1.508257E-04, 2.194942E-03, 2.407634E-02,
2.045719E-02, 3.154309E-03, 1.379132E-04, 2.237403E-06,
2.955046E-07, 2.160090E-05, 4.448143E-04, 1.190621E-02,
2.609831E-02, 1.035142E-02, 1.148483E-03, 4.658676E-05},
        correct  { 4.692559E-02, 6.211871E-04, 3.280219E-07, 3.860131E-15,
0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
8.638024E-03, 4.710140E-02, 1.200938E-02, 8.225027E-06,
4.525178E-10, 4.238455E-15, 0.000000E+00, 0.000000E+00,
4.786808E-04, 1.685266E-02, 4.005842E-02, 2.322158E-03,
5.998091E-06, 2.147037E-09, 1.870239E-13, 5.234564E-18,
3.783126E-06, 3.464375E-04, 5.375396E-03, 2.431827E-02,
2.961144E-03, 3.979691E-05, 1.067401E-07, 8.587682E-11,
3.470559E-08, 4.664788E-06, 1.541758E-04, 5.592549E-03,
6.437990E-03, 8.226465E-04, 1.957708E-05, 1.305188E-07,
2.352498E-10, 5.041676E-08, 2.799380E-06, 3.950756E-04,
2.337875E-03, 1.685040E-03, 2.186062E-04, 7.520523E-06,
9.277760E-13, 3.439243E-10, 3.527198E-08, 1.692636E-05,
3.958404E-04, 1.222261E-03, 6.257662E-04, 7.791164E-05,
2.532257E-15, 1.226521E-12, 2.512334E-10, 4.155520E-07,
2.567706E-05, 2.275790E-04, 3.666844E-04, 1.475488E-04 };
        discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
          //REQUIRE(ranges::equal(sab,correct,equal));
        } // THEN
       // std::cout << (sab | ranges::view::all ) << std::endl;
      } // AND WHEN
      */

    } // WHEN
  } // GIVEN
} // TEST CASE 

























TEST_CASE( "Discrete oscillator treatment (old tests)" ){
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
    std::vector<double> oscEnergies { 0.205, 0.48 }, 
                        oscWeights  { 0.166667, 0.333333 };
    auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);

    WHEN( "alpha and beta valeus are moderately sized" ){
      std::vector<double> 
        alpha = {0.01, 0.04, 0.08, 0.10, 0.40, 0.80, 1.00, 4.00, 8.00, 10.00},
        beta  = {0.00, 1.00, 2.00, 3.00, 4.00},
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

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, oscEnergiesWeights, effectiveTemp, sab );

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sab,sabCorrect,equal));
      } // THEN
    } // WHEN

    WHEN( "alpha and beta both span a huge range" ){
      std::vector<double> 
        alpha = {0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10, 50, 100},
        beta  = {0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10, 50, 100},
        sab {3.783747E+01, 3.389591E+01, 2.424696E+01, 5.632040E-04, 
        5.944319E-05, 2.931679E-05, 3.604359E-05, 1.302505E-05, 2.852185E-10, 
        0.000000E+00, 0.000000E+00, 1.690484E+01, 1.652874E+01, 1.552734E+01, 
        1.827156E+00, 2.486844E-03, 1.482647E-04, 1.794480E-04, 6.472644E-05, 
        7.025305E-09, 0.000000E+00, 0.000000E+00, 1.193902E+01, 1.182827E+01, 
        1.145869E+01, 3.965110E+00, 1.397768E-01, 2.990279E-04, 3.576237E-04, 
        1.289955E-04, 2.781304E-08, 0.000000E+00, 0.000000E+00, 5.288890E+00, 
        5.277051E+00, 5.265238E+00, 4.322478E+00, 2.256523E+00, 1.546854E-03, 
        1.755870E-03, 6.391400E-04, 6.727212E-07, 0.000000E+00, 0.000000E+00, 
        3.696938E+00, 3.693800E+00, 3.690662E+00, 3.379998E+00, 2.474803E+00, 
        3.230450E-03, 3.452499E-03, 1.274852E-03, 2.958873E-06, 5.010057E-30, 
        0.000000E+00, 1.515354E+00, 1.516883E+00, 1.518416E+00, 1.515707E+00, 
        1.457246E+00, 2.191122E-01, 1.593964E-02, 6.432143E-03, 8.011742E-05, 
        1.816200E-23, 0.000000E+00, 9.667891E-01, 9.682187E-01, 9.696371E-01, 
        9.799659E-01, 9.716895E-01, 4.202779E-01, 4.580241E-02, 1.304702E-02, 
        3.455963E-04, 1.457970E-20, 0.000000E+00, 2.128595E-01, 2.133201E-01, 
        2.137239E-01, 2.171325E-01, 2.218032E-01, 2.275933E-01, 1.776032E-01, 
        5.925186E-02, 1.045667E-02, 2.044901E-13, 0.000000E+00, 7.265244E-02, 
        7.282368E-02, 7.299636E-02, 7.441957E-02, 7.611087E-02, 8.707374E-02, 
        9.225780E-02, 7.963887E-02, 3.473308E-02, 3.651351E-10, 1.467265E-24, 
        4.117754E-04, 4.127892E-04, 4.138054E-04, 4.220265E-04, 4.325349E-04, 
        5.257059E-04, 6.655505E-04, 3.201093E-03, 1.151002E-02, 1.287332E-03, 
        6.291107E-10, 1.706488E-06, 1.710763E-06, 1.715049E-06, 1.749688E-06, 
        1.793778E-06, 2.187137E-06, 2.789660E-06, 1.680559E-05, 1.049101E-04, 
        1.531873E-02, 6.550072E-05}, 
        sabCorrect {3.783601E+01, 3.389461E+01, 2.424603E+01, 5.631824E-04, 
        5.944090E-05, 2.931567E-05, 3.604220E-05, 1.302457E-05, 8.703740E-10, 
        0.000000E+00, 0.000000E+00, 1.690159E+01, 1.652556E+01, 1.552435E+01, 
        1.826804E+00, 2.486366E-03, 1.482362E-04, 1.794134E-04, 6.471452E-05, 
        2.158185E-08, 0.000000E+00, 0.000000E+00, 1.193442E+01, 1.182372E+01, 
        1.145428E+01, 3.963585E+00, 1.397231E-01, 2.989129E-04, 3.574861E-04, 
        1.289481E-04, 8.581651E-08, 0.000000E+00, 0.000000E+00, 5.278728E+00, 
        5.266912E+00, 5.255121E+00, 4.314173E+00, 2.252187E+00, 1.543882E-03, 
        1.752496E-03, 6.379641E-04, 2.096558E-06, 1.025453E-18, 0.000000E+00, 
        3.682745E+00, 3.679619E+00, 3.676493E+00, 3.367022E+00, 2.465302E+00, 
        3.218048E-03, 3.439245E-03, 1.270163E-03, 8.560432E-06, 1.821032E-13, 
        0.000000E+00, 1.486489E+00, 1.487989E+00, 1.489492E+00, 1.486835E+00, 
        1.429487E+00, 2.149384E-01, 1.563602E-02, 6.314523E-03, 2.089396E-04, 
        2.440142E-10, 0.000000E+00, 9.303083E-01, 9.316840E-01, 9.330489E-01, 
        9.429879E-01, 9.350238E-01, 4.044191E-01, 4.407416E-02, 1.257785E-02, 
        1.008290E-03, 5.611201E-08, 0.000000E+00, 1.756196E-01, 1.759995E-01, 
        1.763327E-01, 1.791450E-01, 1.829986E-01, 1.877760E-01, 1.465330E-01, 
        4.930239E-02, 2.029639E-02, 3.243924E-05, 2.539787E-12, 4.945904E-02, 
        4.957562E-02, 4.969317E-02, 5.066204E-02, 5.181344E-02, 5.927721E-02, 
        6.280896E-02, 5.479071E-02, 3.622037E-02, 3.282077E-04, 7.837533E-08, 
        6.087852E-05, 6.102842E-05, 6.117870E-05, 6.239442E-05, 6.394848E-05, 
        7.773404E-05, 9.845356E-05, 4.792199E-04, 1.831175E-03, 3.620962E-03, 
        3.530294E-04, 3.795997E-08, 3.805508E-08, 3.815042E-08, 3.892100E-08, 
        3.990205E-08, 4.865977E-08, 6.209764E-08, 3.764874E-07, 2.441230E-06, 
        8.235741E-04, 3.394048E-03};

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, oscEnergiesWeights, effectiveTemp, sab );

      THEN( "scattering law matrix is correctly changed" ){
        REQUIRE(ranges::equal(sab,sabCorrect,equal));
      } // THEN


    } // WHEN
  } // GIVEN

} // TEST CASE
