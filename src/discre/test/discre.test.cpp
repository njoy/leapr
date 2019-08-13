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


  GIVEN( "simplifid water input" ){
    double sc = 0.99187300471670414, scaling = 0.99187300471670414,
    lambda_s = 0.23520419, temp = 296.0,
    twt = 5.5556e-2, tbeta = 0.444444, effectiveTemp = 541.87556285322705;
    std::vector<double> 
      alpha = {0.10, 0.04, 0.08, 0.10, 0.40, 0.80, 1.00, 4.00, 8.00, 10.0},
      beta  = {0.00, 1.00, 2.00, 3.00, 4.00},
      oscEnergies { 0.205,    0.48     },
      oscWeights  { 0.166667, 0.666666 },
      sab { 11.98804915, 3.610205E-4, 4.864924E-4, 4.535244E-4, 1.156846E-4, 5.951624125, 1.451989E-3, 1.926084E-3, 1.791418E-3, 4.772776E-4, 4.169703384, 2.909939E-3, 3.809574E-3, 3.539618E-3, 9.786056E-4, 3.712543196, 3.638115E-3, 4.737517E-3, 4.401010E-3, 1.235316E-3, 1.739126443, 1.435423E-2, 1.771050E-2, 1.652928E-2, 5.379703E-3, 1.133534359, 3.434249E-2, 3.264296E-2, 3.082858E-2, 1.137667E-2, 0.974899977, 5.136561E-2, 3.922902E-2, 3.728513E-2, 1.447020E-2, 0.289595362, 0.203273979, 9.854002E-2, 9.147211E-2, 5.420127E-2, 0.114055698, 0.133646508, 0.109956581, 9.825605E-2, 7.232408E-2, 7.821841E-2, 0.101069704, 9.734660E-2, 8.954631E-2, 7.085870E-2 },
      correctSAB { 11.9834763629, 3.60882816E-4, 4.86306841E-4, 4.53351406E-4, 1.15640497E-4, 5.94254842792, 1.44977504E-3, 1.92314734E-3, 1.78868636E-3, 4.76549837E-4, 4.15699622913, 2.90107127E-3, 3.79796452E-3, 3.52883189E-3, 9.75625188E-4, 3.69840613633, 3.62426163E-3, 4.71947768E-3, 4.38425293E-3, 1.23064037E-3, 1.71278759424, 1.41368548E-2, 1.74423270E-2, 1.62791889E-2, 5.29930548E-3, 1.09945995083, 3.33101949E-2, 3.16618467E-2, 2.99023380E-2, 1.10362578E-2, 0.93840644484, 4.94428903E-2, 3.77607212E-2, 3.58899321E-2, 1.39300603E-2, 0.24860853768, 0.17450439559, 8.45937177E-2, 7.85264800E-2, 4.65315819E-2, 8.40558237E-2, 9.84937256E-2, 8.10350275E-2, 7.24123401E-2, 5.33019471E-2, 5.34101406E-2, 6.90137869E-2, 6.64716132E-2, 6.11455485E-2, 4.83856381E-2 };


    //auto oscEnergies = ranges::view::iota(1,5); 
    //auto oscWeights  = ranges::view::iota(2,6); 
    auto oscEnergiesWeights = ranges::view::zip(oscEnergies,oscWeights);
//    std::cout << std::endl;

    discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, temp, oscEnergiesWeights, effectiveTemp, sab );
 //   std::cout << std::endl;
      THEN( "scattering law matrix is correctly changed" ){
        //REQUIRE(ranges::equal(sab,correctSAB,equal));
      } // THEN

    //std::cout << (sab | ranges::view::all) << std::endl;
    /*
      std::vector<double> osc_energies{0.205, 0.48}, osc_weights{0.166667, 0.666666},
      alpha = {0.10, 0.04, 0.08, 0.10, 0.40, 0.80, 1.00, 4.00, 8.00, 10.0},
      beta  = {0.00, 1.00, 2.00, 3.00, 4.00};
      double t_eff = 81178.935219;

      //std::vector<double> sym_sab (alpha.size()*beta.size(),0.0);
      //for ( size_t i = 0; i < sym_sab.size(); ++i ){ sym_sab[i] = i+1; }

      std::vector<double> sym_sab  =  { 11.98804915, 3.610205E-4, 4.864924E-4, 4.535244E-4, 1.156846E-4, 5.951624125, 1.451989E-3, 1.926084E-3, 1.791418E-3, 4.772776E-4, 4.169703384, 2.909939E-3, 3.809574E-3, 3.539618E-3, 9.786056E-4, 3.712543196, 3.638115E-3, 4.737517E-3, 4.401010E-3, 1.235316E-3, 1.739126443, 1.435423E-2, 1.771050E-2, 1.652928E-2, 5.379703E-3, 1.133534359, 3.434249E-2, 3.264296E-2, 3.082858E-2, 1.137667E-2, 0.974899977, 5.136561E-2, 3.922902E-2, 3.728513E-2, 1.447020E-2, 0.289595362, 0.203273979, 9.854002E-2, 9.147211E-2, 5.420127E-2, 0.114055698, 0.133646508, 0.109956581, 9.825605E-2, 7.232408E-2, 7.821841E-2, 0.101069704, 9.734660E-2, 8.954631E-2, 7.085870E-2 };

      std::cout << sym_sab.size() << std::endl;
      double lambda_s = 2.2941534E-3, sc = 1, scaling = 1, tbeta = 2, twt = 0.3;

    auto oscEnergiesWeights = ranges::view::zip(osc_energies,osc_weights);

      //discre( sc, scaling, lambda_s, twt, tbeta, alpha, beta, 
       //       temp, oscEnergiesWeights, t_eff, sym_sab );

      std::vector<double> correctSymSab {0.9575582, 1.914953, 2.872356, 
        3.829659, 4.807700, 5.501617, 6.418366, 7.335187, 8.251386, 9.253313, 
        9.256188, 10.09717, 10.93859, 11.77685, 12.85839, 11.36789, 12.07716, 
        12.78842, 13.48618, 14.74575, 10.74546, 11.25461, 11.77105, 12.24116, 
        13.75305};

      THEN( "scattering law matrix is correctly changed" ){
        //REQUIRE(ranges::equal(sym_sab,correctSymSab,equal));
      } // THEN

      */

  } // GIVEN


} // TEST CASE













