#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "leapr.cpp"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

TEST_CASE( "leapr" ){
  int nphon, ncold, lat;
  double sps, awr, aws, delta, twt, c, tbeta, dka, cfrac;
  std::vector<double> alpha, beta, temp, rho, oscE, oscW, kappa;

  GIVEN( "coarse alpha, beta grids (for quick testing)" ) {
    WHEN( "continuous, translational, and discrete oscillator options used" ) {
      nphon = 100;
      awr   = 0.99917;  ncold = 0; 
      aws   = 15.85316; sps = 3.8883; 
      lat   = 0;
      alpha = {.01008,   .033,     0.050406, 0.504060, 0.554466, 1.873790, 
               2.055660, 5.671180, 7.472890, 20.52030, 32.46250, 50.0};
      beta  = {0.000000, 0.006375, 0.038250, 0.06575, 0.564547, 1.77429,
               3.422470, 8.258620, 18.72180, 25.0};
      temp  = { 296.0 };                                          
      delta = 0.00255;                                           
      rho = {0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
             0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
             0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
             0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 
             0.0839, 0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 
             0.06652, 0.065, 0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 
             0.05462, 0.0535, 0.0525, 0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 
             0.0459, 0.04478, 0.04366, 0.04288, 0.04244, 0.042, 0.0, };
      twt    = 0.055556;   c = 0.0;   tbeta = 0.444444;
      oscE   = { 0.205,    0.48};                                 
      oscW   = { 0.166667, 0.333333 };                            
      dka    = 0.0;
      cfrac  = 0.0;
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);

      /*
      AND_WHEN ( "alpha and beta scaling not needed (lat = 0)" ){
        lat   = 0;            
        auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac );
  
        std::vector<double> sabCorrect { 1.188670E+1, 1.168376E+1, 6.300891E+0, 
        1.781912E+0, 3.162210E-4, 4.692060E-4, 3.081720E-4, 2.514578E-5, 
        1.701154E-4, 6.09094E-10, 6.527907E+0, 6.500191E+0, 5.437262E+0, 
        3.734029E+0, 1.050239E-3, 1.518856E-3, 9.981392E-4, 1.995450E-4, 
        1.200483E-3, 1.282587E-8, 5.256809E+0, 5.241938E+0, 4.691417E+0, 
        3.687751E+0, 1.615068E-3, 2.303674E-3, 1.516466E-3, 3.464544E-4, 
        1.860820E-3, 3.807045E-8, 1.478889E+0, 1.480819E+0, 1.486305E+0, 
        1.470555E+0, 1.288961E-1, 2.017464E-2, 1.408493E-2, 7.462729E-3, 
        1.027139E-2, 1.373864E-5, 1.392511E+0, 1.394460E+0, 1.401735E+0, 
        1.388500E+0, 1.549325E-1, 2.190231E-2, 1.538571E-2, 8.382694E-3, 
        1.082756E-2, 1.748288E-5, 5.552458E-1, 5.564724E-1, 5.626823E-1, 
        5.679902E-1, 3.628562E-1, 5.450169E-2, 4.333900E-2, 2.265130E-2, 
        1.672152E-2, 3.550178E-4, 5.089281E-1, 5.100789E-1, 5.158600E-1, 
        5.210213E-1, 3.557392E-1, 5.783144E-2, 4.633066E-2, 2.381095E-2, 
        1.691639E-2, 4.420672E-4, 1.455618E-1, 1.459800E-1, 1.481698E-1, 
        1.499734E-1, 1.599368E-1, 9.192485E-2, 7.299081E-2, 3.854093E-2, 
        1.533922E-2, 3.925552E-3, 9.083275E-2, 9.109755E-2, 9.245460E-2, 
        9.367509E-2, 1.060097E-1, 8.399707E-2, 7.119046E-2, 4.314996E-2, 
        1.459694E-2, 6.322492E-3, 6.735360E-3, 6.756551E-3, 6.863877E-3, 
        6.958397E-3, 8.741588E-3, 1.328617E-2, 1.925801E-2, 3.149450E-2, 
        2.194943E-2, 1.818368E-2, 9.378210E-4, 9.408032E-4, 9.559292E-4, 
        9.692905E-4, 1.230320E-3, 2.054059E-3, 3.601626E-3, 1.031145E-2, 
        1.933435E-2, 1.935657E-2, 6.228358E-5, 6.248242E-5, 6.348621E-5, 
        6.436519E-5, 8.208659E-5, 1.425013E-4, 2.777029E-4, 1.245381E-3, 
        6.606240E-3, 1.044380E-2};
  
        double debyeWaller   = 0.2352041964;
        std::vector<double> effectiveTemps {1397.671178314};

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( debyeWaller == Approx(std::get<2>(out)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
      */
      AND_WHEN ( "alpha and beta scaling is requested (lat = 1)" ){
        lat   = 1;            
        auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac );
        std::vector<double> sabCorrect { 1.193557E+1, 1.173342E+1, 6.361294E+0, 
          1.816275E+0, 3.124950E-4, 4.639237E-4, 3.122224E-4, 7.882191E-5, 
          1.011487E-5,7.888274E-10, 6.555078E+0, 6.527275E+0, 5.467981E+0, 
          3.766171E+0, 1.038090E-3, 1.501995E-3, 1.011017E-3, 5.830163E-4, 
          8.176285E-5, 1.597145E-8, 5.278888E+0, 5.263947E+0, 4.715237E+0, 
          3.714425E+0, 1.596591E-3, 2.278303E-3, 1.535826E-3, 9.525713E-4, 
          1.452123E-4, 4.673899E-8, 1.486402E+0, 1.488315E+0, 1.493827E+0, 
          1.477974E+0, 1.316516E-1, 1.997813E-2, 1.423771E-2, 1.018858E-2, 
          4.169200E-3, 1.556686E-5, 1.399711E+0, 1.401644E+0, 1.408946E+0, 
          1.395614E+0, 1.581786E-1, 2.169138E-2, 1.555017E-2, 1.108814E-2, 
          4.735061E-3, 1.974070E-5, 5.593066E-1, 5.605274E-1, 5.667146E-1, 
          5.720510E-1, 3.663735E-1, 5.415787E-2, 4.367399E-2, 2.386446E-2, 
          1.252101E-2, 3.835034E-4, 5.127808E-1, 5.139310E-1, 5.196984E-1, 
          5.248412E-1, 3.591896E-1, 5.751228E-2, 4.667544E-2, 2.488494E-2, 
          1.293548E-2, 4.758096E-4, 1.473396E-1, 1.477655E-1, 1.499618E-1, 
          1.517737E-1, 1.616866E-1, 9.247084E-2, 7.336139E-2, 3.867818E-2, 
          1.372647E-2, 4.048636E-3, 9.213506E-2, 9.240120E-2, 9.376556E-2, 
          9.499839E-2, 1.073307E-1, 8.467699E-2, 7.158413E-2, 4.326283E-2, 
          1.351358E-2, 6.455035E-3, 6.930727E-3, 6.952350E-3, 7.061863E-3, 
          7.158157E-3, 8.973773E-3, 1.358021E-2, 1.959090E-2, 3.179135E-2, 
          2.203982E-2, 1.822503E-2, 9.778850E-4, 9.809690E-4, 9.966260E-4, 
          1.010448E-3, 1.279979E-3, 2.127476E-3, 3.710551E-3, 1.050751E-2, 
          1.954142E-2, 1.944808E-2, 6.624281E-5, 6.645256E-5, 6.751134E-5, 
          6.843838E-5, 8.710059E-5, 1.505165E-4, 2.916502E-4, 1.290367E-3, 
          6.740754E-3, 1.056903E-2};
        double debyeWaller   = 0.2352041964;
        std::vector<double> effectiveTemps {1397.671178314};

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( debyeWaller == Approx(std::get<2>(out)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
      AND_WHEN( "material is cold" ){
        temp  = { 100.0 };                                          
        auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac );

        std::vector<double> sabCorrect { 1.191059E+1, 1.170723E+1, 6.313288E+0, 
          1.785018E+0, 1.591388E-5, 1.977783E-5, 3.162744E-5, 5.709636E-5, 
          8.766729E-6, 2.004754E-9, 6.570177E+0, 6.542268E+0, 5.472154E+0, 
          3.757436E+0, 5.258778E-5, 6.470506E-5, 1.031839E-4, 1.860698E-4, 
          2.855014E-5, 3.566762E-8, 5.308403E+0, 5.293371E+0, 4.737147E+0, 
          3.723126E+0, 8.070512E-5, 9.875910E-5, 1.572592E-4, 2.834417E-4, 
          4.352879E-5, 8.395187E-8, 1.616967E+0, 1.619044E+0, 1.624822E+0, 
          1.607225E+0, 1.252489E-1, 9.603589E-4, 1.504800E-3, 2.708091E-3, 
          4.368001E-4, 3.862403E-5, 1.535347E+0, 1.537456E+0, 1.545270E+0, 
          1.530291E+0, 1.539718E-1, 1.052854E-3, 1.647889E-3, 2.966128E-3, 
          4.811426E-4, 4.881106E-5, 7.499578E-1, 7.516101E-1, 7.599277E-1, 
          7.668914E-1, 4.630478E-1, 4.192283E-3, 4.983860E-3, 9.056953E-3, 
          1.687982E-3, 5.798369E-4, 7.055258E-1, 7.071121E-1, 7.150962E-1, 
          7.220613E-1, 4.659443E-1, 5.252341E-3, 5.387029E-3, 9.805895E-3, 
          1.860422E-3, 7.168131E-4, 3.175353E-1, 3.183719E-1, 3.225950E-1, 
          3.262726E-1, 3.279453E-1, 6.996016E-2, 1.131200E-2, 2.105390E-2, 
          5.444883E-3, 6.518448E-3, 2.396012E-1, 2.402501E-1, 2.435240E-1, 
          2.463885E-1, 2.625363E-1, 9.412094E-2, 1.390007E-2, 2.456444E-2, 
          7.245434E-3, 9.318142E-3, 5.184109E-2, 5.199249E-2, 5.275752E-2, 
          5.342847E-2, 6.426167E-2, 6.634176E-2, 3.319597E-2, 2.811754E-2, 
          1.693344E-2, 1.658023E-2, 1.635642E-2, 1.640538E-2, 1.665411E-2, 
          1.687429E-2, 2.082778E-2, 2.691698E-2, 2.342642E-2, 1.982073E-2, 
          1.878867E-2, 1.665706E-2, 3.443731E-3, 3.454243E-3, 3.507459E-3, 
          3.554265E-3, 4.448901E-3, 6.591632E-3, 8.198762E-3, 9.342037E-3, 
          1.423867E-2, 1.370927E-2};
        double debyeWaller   = 5.68713599E-2;
        std::vector<double> effectiveTemps {1366.18152086};

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( debyeWaller == Approx(std::get<2>(out)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
      AND_WHEN( "material is hot" ){
        temp  = { 600.0 };                                          
        auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac );

        std::vector<double> sabCorrect { 1.182050E+1, 1.161877E+1, 6.267602E+0, 
          1.775094E+0, 2.538069E-3, 1.274758E-3, 3.453815E-6, 2.589951E-7, 
          2.71834E-15, 0.000000E+0, 6.414739E+0, 6.387587E+0, 5.344998E+0, 
          3.674163E+0, 8.254109E-3, 4.092687E-3, 3.690549E-5, 2.840638E-6, 
          5.29688E-13, 3.88683E-20, 5.121307E+0, 5.106933E+0, 4.572468E+0, 
          3.597790E+0, 1.253030E-2, 6.194625E-3, 8.554847E-5, 6.672075E-6, 
          2.151270E-6, 1.84458E-13, 1.207943E+0, 1.209742E+0, 1.214528E+0, 
          1.204356E+0, 1.856763E-1, 5.413021E-2, 8.248944E-3, 8.973797E-4, 
          1.418359E-4, 8.605015E-8, 1.119380E+0, 1.120947E+0, 1.126974E+0, 
          1.119192E+0, 2.084083E-1, 5.868859E-2, 9.891214E-3, 1.101023E-3, 
          1.650710E-4, 1.300575E-7, 3.296997E-1, 3.309046E-1, 3.366573E-1, 
          3.392579E-1, 2.810407E-1, 1.296941E-1, 5.873684E-2, 8.093005E-3, 
          8.179071E-4, 2.181782E-5, 2.937824E-1, 2.948494E-1, 3.003364E-1, 
          3.026732E-1, 2.655510E-1, 1.336042E-1, 6.440555E-2, 9.156526E-3, 
          9.122884E-4, 3.149881E-5, 6.493049E-2, 6.513407E-2, 6.616440E-2, 
          6.709328E-2, 8.129466E-2, 9.848076E-2, 9.553444E-2, 3.456580E-2, 
          4.882244E-3, 1.157556E-3, 3.851669E-2, 3.863880E-2, 3.925635E-2, 
          3.976768E-2, 4.933831E-2, 6.852243E-2, 7.992601E-2, 4.483372E-2, 
          8.770187E-3, 2.590310E-3, 1.913698E-3, 1.919809E-3, 1.950662E-3, 
          1.977712E-3, 2.520742E-3, 4.266602E-3, 7.722683E-3, 2.172338E-2, 
          2.834815E-2, 2.092619E-2, 1.642067E-4, 1.647312E-4, 1.673778E-4, 
          1.696941E-4, 2.169625E-4, 3.781273E-4, 7.515546E-4, 3.461227E-3, 
          1.537487E-2, 2.126633E-2, 5.110181E-6, 5.126508E-6, 5.208941E-6, 
          5.281167E-6, 6.762964E-6, 1.198200E-5, 2.513165E-5, 1.576405E-4, 
          2.127419E-3, 5.732186E-3 };
        double debyeWaller = 0.78234716272468541;
        std::vector<double> effectiveTemps {1507.2963449849456};
          
        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( debyeWaller == Approx(std::get<2>(out)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN

   /*
  GIVEN( "test9 (simple H in H2O)" ) {
    WHEN( "continuous, translational, and discrete oscillator options used" ) {
      nphon = 100;
      awr   = 0.99917;  ncold = 0; 
      aws   = 15.85316; sps = 3.8883; 
      lat   = 1;
      alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 0.100812, 
      0.151218, 0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 
      0.504060, 0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 
      0.976435, 1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 
      1.873790, 2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 
      3.578320, 3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 
      6.816500, 7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 
      12.96740, 14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 
      24.65260, 27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 
      46.85030, 50.0 };
      beta = { 0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 
      0.065750, 0.0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 
      0.483897, 0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 
      1.048440, 1.129090, 1.209740, 1.290390, 1.371040, 1.451690, 1.532340, 
      1.612990, 1.693640, 1.774290, 1.854940, 1.935590, 2.016240, 2.096890, 
      2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 2.580790, 2.669500, 
      2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 3.595360, 
      3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 
      5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 
      9.637220, 10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 
      17.17330, 18.72180, 20.42450, 22.29760, 24.35720, 25.0 };


      temp  = { 296.0 };                                          
      delta = 0.00255;                                           
      rho = {0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
             0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
             0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
             0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 
             0.0839, 0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 
             0.06652, 0.065, 0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 
             0.05462, 0.0535, 0.0525, 0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 
             0.0459, 0.04478, 0.04366, 0.04288, 0.04244, 0.042, 0.0, };
      twt    = 0.055556;   c = 0.0;   tbeta = 0.444444;
      oscE   = { 0.205,    0.48};                                 
      oscW   = { 0.166667, 0.333333 };                            
      dka    = 0.0;
      cfrac  = 0.0;
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);


      auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                        rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac );
      std::vector<double> sab = std::get<0>(out);
      std::vector<double> sab_0_99 { 1.193557E+1, 1.173342E+1, 1.115293E+1, 9.042512E+0, 6.361295E+0, 3.862215E+0, 1.816276E+0, 6.976093E-1, 1.990600E-2, 5.694252E-4, 3.233509E-4, 2.982178E-4, 2.996330E-4, 2.961039E-4, 3.053551E-4, 3.177549E-4, 3.265274E-4, 3.366876E-4, 3.482198E-4, 3.577151E-4, 3.690864E-4, 3.825873E-4, 3.943943E-4, 4.067551E-4, 4.197706E-4, 4.318199E-4, 4.422418E-4, 4.513514E-4, 4.592862E-4, 4.676470E-4, 4.766943E-4, 4.858306E-4, 4.954565E-4, 5.093038E-4, 5.232462E-4, 5.360422E-4, 5.478218E-4, 5.586641E-4, 5.670666E-4, 5.670803E-4, 5.574108E-4, 5.388335E-4, 5.099740E-4, 4.642119E-4, 4.116225E-4, 3.618202E-4, 3.183550E-4, 2.801315E-4, 2.468202E-4, 2.184691E-4, 1.934883E-4, 1.714009E-4, 1.502529E-4, 1.302423E-4, 1.135101E-4, 9.827482E-5, 8.024491E-5, 1.829132E-5, 2.490368E-7, 1.936146E-7, 3.065288E-7, 1.477736E-7, 1.342006E-7, 1.349326E-7, 8.946933E-8, 4.469996E-8, 2.540386E-8, 8.788727E-9, 2.848220E-9, 1.400958E-8, 4.429853E-8, 7.649074E-8, 6.102628E-8, 2.022216E-8, 1.531153E-8, 9.770919E+0, 9.675363E+0, 9.353680E+0, 8.135114E+0, 6.428072E+0, 4.614642E+0, 2.783247E+0, 1.467645E+0, 1.338837E-1, 5.248679E-3, 4.879553E-4, 4.447880E-4, 4.446880E-4, 4.413456E-4, 4.539861E-4, 4.719887E-4, 4.855082E-4, 5.005741E-4, 5.173926E-4, 5.318829E-4, 5.487859E-4, 5.685026E-4, 5.863043E-4, 6.046558E-4, 6.238830E-4 };

      checkPartOfVec( sab, sab_0_99 );


    } // WHEN
 } // GIVEN
  */
} // TEST CASE


