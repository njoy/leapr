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
        std::vector<double> 
          alpha = {0.1, 1.0, 10},
          beta  = {0.0, 5.0, 50},
          sab { 6.28055E-3, 1.48422E-3, 1.78104E-29, 5.10424E-2, 1.52272E-2, 
          4.7449E-20, 7.01340E-2, 8.69301E-2, 8.7082E-10}, 
          correct { 6.25390E-3, 1.47793E-3, 0.00000000, 4.89179E-2, 1.45936E-2, 
          3.1085E-11, 4.58457E-2, 5.68386E-2, 8.52874E-5 };
          discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                  oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
           REQUIRE(ranges::equal(sab,correct,equal));
        } // THEN
      } // AND WHEN
    }
    /*
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
        correct { 6.42751E-7, 6.42911E-7, 6.42780E-7, 6.42911E-7, 6.43074E-7, 
        6.44381E-7, 6.46015E-7, 6.59084E-7, 3.21332E-5, 3.21412E-5, 3.21347E-5, 
        3.21412E-5, 3.21494E-5, 3.22147E-5, 3.22964E-5, 3.29497E-5, 6.42735E-6, 
        6.42895E-6, 6.42764E-6, 6.42895E-6, 6.43058E-6, 6.44365E-6, 6.45999E-6, 
        6.59067E-6, 3.21332E-5, 3.21412E-5, 3.21347E-5, 3.21412E-5, 3.21494E-5, 
        3.22147E-5, 3.22964E-5, 3.29497E-5, 6.42576E-5, 6.42736E-5, 6.42606E-5, 
        6.42736E-5, 6.42900E-5, 6.44206E-5, 6.45839E-5, 6.58905E-5, 3.20936E-4, 
        3.21016E-4, 3.20951E-4, 3.21016E-4, 3.21098E-4, 3.21750E-4, 3.22566E-4, 
        3.29091E-4, 6.40995E-4, 6.41154E-4, 6.41024E-4, 6.41154E-4, 6.41317E-4, 
        6.42620E-4, 6.44249E-4, 6.57280E-4, 3.17006E-3, 3.17085E-3, 3.17020E-3, 
        3.17085E-3, 3.17165E-3, 3.17809E-3, 3.18614E-3, 3.25055E-3 };
        discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                oscEnergiesWeights, effectiveTemp, sab );
        THEN( "scattering law matrix is correctly changed" ){
          REQUIRE(ranges::equal(sab,correct,equal));
        } // THEN
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
      AND_WHEN( "Rather large alpha, beta values" ){
        std::vector<double> 
        alpha1 = {60, 100, 120},
        beta1  = {0,  100, 120},
        sab1 { 1.70896E-4, 5.21119E-8, 7.0407E-11, 2.39298E-6, 1.37913E-4, 
        2.23740E-6, 2.95504E-7, 1.14848E-3, 4.65867E-5 },
        correct1 { 1.36280E-8, 1.17072E-7, 2.12854E-8, 1.1482E-13, 2.15253E-6, 
        4.70458E-6, 8.8005E-16, 1.16734E-6, 5.99689E-6 },

        alpha2 {60, 100, 120, 150},
        beta2  {0,  100, 120, 150},
        sab2 { 1.70896E-4, 5.21119E-8, 7.0407E-11, 8.5425E-16, 2.39298E-6, 
        1.37913E-4, 2.23740E-6, 1.04788E-9, 2.95504E-7, 1.14848E-3, 4.65867E-5, 
        8.40991E-8, 1.31808E-8, 7.62407E-3, 9.66387E-4, 9.37698E-6},
        correct2 { 1.36280E-8, 1.17072E-7, 2.12854E-8, 1.1333E-11, 1.1482E-13, 
        2.15253E-6, 4.70458E-6, 4.31086E-7, 8.8005E-16, 1.16734E-6, 5.99689E-6, 
        4.54521E-6, 8.7052E-19, 1.97621E-7, 2.80726E-6, 2.03141E-5},

        alpha3 {10,30,60,100,120,150},
        beta3  {0, 50,100,120,150},
        sab3 { 7.01340E-2, 8.7082E-10, 1.0288E-23, 2.3337E-30, 0.00000000, 
        5.38240E-3, 4.48589E-5, 8.8775E-14, 5.0053E-18, 5.6793E-25, 1.70896E-4, 
        6.77023E-3, 5.21119E-8, 7.0407E-11, 8.5425E-16, 2.39298E-6, 2.85699E-2, 
        1.37913E-4, 2.23740E-6, 1.04788E-9, 2.95504E-7, 2.26650E-2, 1.14848E-3, 
        4.65867E-5, 8.40991E-8, 1.31808E-8, 8.35859E-3, 7.62407E-3, 9.66387E-4, 
        9.37698E-6},
        correct3 { 6.17897E-3, 1.34999E-8, 0.00000000, 0.000000000, 0.00000000, 
        2.36273E-5, 2.10282E-4, 3.7982E-10, 1.2303E-13, 0.00000000, 1.74719E-8, 
        4.57293E-4, 8.64524E-6, 7.73706E-8, 1.1333E-11, 1.8790E-13, 8.18326E-6, 
        3.00123E-4, 3.55855E-5, 4.32883E-7, 1.0024E-15, 3.77351E-7, 2.03275E-4, 
        6.32648E-5, 4.70616E-6, 8.7052E-19, 2.06917E-9, 4.03784E-5, 4.21806E-5, 
        2.22134E-5};

        discre( sc, scaling, lambda_s, twt, tbeta, alpha1, beta1, temp, 
                oscEnergiesWeights, effectiveTemp, sab1 );
        discre( sc, scaling, lambda_s, twt, tbeta, alpha2, beta2, temp, 
                oscEnergiesWeights, effectiveTemp, sab2 );
        discre( sc, scaling, lambda_s, twt, tbeta, alpha3, beta3, temp, 
                oscEnergiesWeights, effectiveTemp, sab3 );

        THEN( "scattering law matrix is correctly changed" ){
          REQUIRE(ranges::equal(sab1,correct1,equal));
          REQUIRE(ranges::equal(sab2,correct2,equal));
          REQUIRE(ranges::equal(sab3,correct3,equal));
        } // THEN
      } // AND WHEN
    } // WHEN
    */
      /*
    WHEN( "Hefty number of oscillators present" ){
      std::vector<double> oscEnergies {0.01, 0.02, 0.05, 0.10, 0.20, 0.50}, 
                          oscWeights  {0.10, 0.10, 0.10, 0.10, 0.05, 0.05};
      auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);
      std::vector<double> 
        alpha {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0},
        beta  {0.0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0},
        sab { 
6.280552E-03, 3.474793E-03, 3.205773E-03, 3.990279E-03,
5.469996E-03, 1.484229E-03, 3.713332E-06, 7.253001E-09,
1.227414E-02, 6.810920E-03, 6.290901E-03, 7.829325E-03,
1.073361E-02, 2.982076E-03, 1.509490E-05, 5.914239E-08,
2.863284E-02, 1.602976E-02, 1.485720E-02, 1.848317E-02,
2.534654E-02, 7.537271E-03, 9.852097E-05, 9.737234E-07,
5.104245E-02, 2.899986E-02, 2.703259E-02, 3.360907E-02,
4.611293E-02, 1.522726E-02, 4.175489E-04, 8.379662E-06,
8.120994E-02, 4.751448E-02, 4.479115E-02, 5.562597E-02,
7.641469E-02, 3.032828E-02, 1.800458E-03, 7.449372E-05,
1.033833E-01, 6.596428E-02, 6.420017E-02, 7.954266E-02,
1.098312E-01, 6.628287E-02, 1.166764E-02, 1.316039E-03,
7.013400E-02, 5.127893E-02, 5.229343E-02, 6.471524E-02,
9.047949E-02, 8.693012E-02, 3.705110E-02, 9.466065E-03,
3.801206E-02, 3.131686E-02, 3.316846E-02, 4.110206E-02,
5.837130E-02, 7.690541E-02, 5.641066E-02, 2.396173E-02},

correct {
5.676546E-03, 3.371587E-03, 1.642089E-01, 3.941329E-03,
5.318661E-03, 1.499506E-03, 8.467173E-06, 9.756671E-08,
1.013705E-02, 6.426192E-03, 2.711905E-01, 7.642262E-03,
1.017326E-02, 3.058140E-03, 3.578610E-05, 5.959609E-07,
1.908606E-02, 1.404602E-02, 3.936910E-01, 1.742441E-02,
2.251829E-02, 8.103712E-03, 1.789460E-04, 1.975495E-06,
2.711094E-02, 2.293365E-02, 3.518622E-01, 2.978069E-02,
3.771458E-02, 1.718466E-02, 7.508540E-04, 1.559130E-05,
3.293573E-02, 3.087408E-02, 1.953727E-01, 4.213642E-02,
5.442590E-02, 3.456911E-02, 3.079567E-03, 1.379333E-04,
2.541736E-02, 2.573976E-02, 5.049130E-02, 3.652988E-02,
5.142962E-02, 5.867942E-02, 1.673464E-02, 2.321456E-03,
6.998583E-03, 7.409937E-03, 8.798037E-03, 1.095410E-02,
1.699866E-02, 2.890906E-02, 2.061739E-02, 7.001931E-03,
1.485627E-03, 1.627021E-03, 1.950335E-03, 2.529318E-03,
4.188025E-03, 1.042605E-02, 1.485453E-02, 9.618843E-03 };
    

      alpha = {0.1, 0.5, 1.0, 5.0, 10.0, 15.0};
      beta  = {0.0, 0.5, 1.0, 5.0, 10.0, 15.0};
      sab = {6.280552E-03, 3.205773E-03, 3.990279E-03, 1.484229E-03,
3.713332E-06, 7.253001E-09, 2.863284E-02, 1.485720E-02,
1.848317E-02, 7.537271E-03, 9.852097E-05, 9.737234E-07,
5.104245E-02, 2.703259E-02, 3.360907E-02, 1.522726E-02,
4.175489E-04, 8.379662E-06, 1.033833E-01, 6.420017E-02,
7.954266E-02, 6.628287E-02, 1.166764E-02, 1.316039E-03,
7.013400E-02, 5.229343E-02, 6.471524E-02, 8.693012E-02,
3.705110E-02, 9.466065E-03, 3.801206E-02, 3.316846E-02,
4.110206E-02, 7.690541E-02, 5.641066E-02, 2.396173E-02};
      correct = {5.709323E-03, 1.320894E-01, 3.909177E-03, 1.476674E-03,
7.188199E-06, 9.596417E-08, 1.943832E-02, 3.188184E-01,
1.708664E-02, 7.608528E-03, 1.774248E-04, 1.975495E-06,
2.755194E-02, 2.880581E-01, 2.926628E-02, 1.554017E-02,
7.503221E-04, 1.559097E-05, 2.429154E-02, 4.575932E-02,
3.575106E-02, 5.079566E-02, 1.668751E-02, 2.302610E-03,
6.647755E-03, 8.497614E-03, 1.098011E-02, 2.608253E-02,
2.043681E-02, 7.004838E-03, 1.409505E-03, 1.873685E-03,
2.454772E-03, 9.714565E-03, 1.463676E-02, 9.617850E-03};


      alpha = {5.0, 10.0, 15.0};
      beta  = {5.0, 10.0, 15.0};
      sab = {6.628287E-02, 1.166764E-02, 1.316039E-03, 8.693012E-02,
3.705110E-02, 9.466065E-03, 7.690541E-02, 5.641066E-02,
2.396173E-02 };
      correct = {4.833928E-02, 1.667136E-02, 2.266506E-03, 3.183675E-02,
2.060186E-02, 7.004407E-03, 1.716947E-02, 1.538098E-02,
9.627496E-03};


      std::cout << std::endl;

      discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, 
                oscEnergiesWeights, effectiveTemp, sab );
      THEN( "scattering law matrix is correctly changed" ){
//        REQUIRE(ranges::equal(sab,correct,equal));
      } // THEN

    std::cout << std::endl;
    std::cout << (sab|ranges::view::all) << std::endl;
    std::cout << std::endl;

    } // WHEN
    */

  } // GIVEN
} // TEST CASE 














/*











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
*/
