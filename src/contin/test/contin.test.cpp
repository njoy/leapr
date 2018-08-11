#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin.h"


void checkSabLambdaTeff( const std::vector<double>& correctSab, 
    const std::tuple<double,double>& output, 
    const Eigen::Tensor<double,3>& sab,  
    const double& lambda, const double& teff, const double& tol ){

  REQUIRE( sab.dimension(0)*sab.dimension(1)*sab.dimension(2) == correctSab.size() );

  int l = 0;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      for ( int k = 0; k < sab.dimension(2); ++k ){
        REQUIRE( sab(i,j,k) == Approx(correctSab[l]).epsilon(tol) );
	l += 1;
      }
    }
  }
  REQUIRE( std::get<0>(output) == Approx(lambda).epsilon(1e-6) );
  REQUIRE( std::get<1>(output) == Approx(teff).epsilon(1e-6) );
}


TEST_CASE( "contin eigen" ){

  int ntempr, nphon, itemp;
  double delta, tbeta, tev, sc, scaling, lambda_s, t_eff;
  std::vector<double> alpha, beta, rho, expected;
  std::tuple<double,double> output;

  GIVEN( "input values from input card and leapr subroutine" ){
    ntempr = 1; nphon = 3; itemp = 0;
    delta = 0.1; tbeta = 1.0; tev = 0.01723477; sc = 0.0253/tev; scaling = sc;
    alpha = { 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28};
    beta  = { 0.00, 0.15, 0.30, 0.60, 1.20 };
    rho   = { 0.002, 0.004, 0.02, 0.04, 0.2, 0.4 };


    WHEN( "3rd order expansion, with alpha & beta vals scaled by 0.0253/tev" ){

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      THEN( "contin output matches expected value" ){
        lambda_s = 4.38473153E-2, t_eff = 12.279863466;
        expected = {6.474331963E-7, 7.658564709E-7, 8.842797456E-7, 
          1.121126294E-6, 1.594819393E-6, 1.294036143E-6, 1.530731302E-6, 
          1.767426462E-6, 2.240816781E-6, 3.187597419E-6, 2.584754483E-6, 
          3.057539823E-6, 3.530325163E-6, 4.475895844E-6, 6.367037205E-6, 
          5.156263261E-6, 6.099418664E-6, 7.042574067E-6, 8.928884872E-6, 
          1.270150648E-5, 1.025974716E-5, 1.213643384E-5, 1.401312052E-5, 
          1.776649387E-5, 2.527324059E-5, 2.030999488E-5, 2.402516585E-5, 
          2.774033682E-5, 3.517067875E-5, 5.003136263E-5, 3.979478204E-5, 
          4.707463949E-5, 5.435449694E-5, 6.891421183E-5, 9.803364162E-5, 
          7.638863877E-5, 9.036457859E-5, 1.043405184E-4, 1.322923980E-4, 
          1.881961573E-4};
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff, 1.0e-6 );
      } // THEN
    } // WHEN

    WHEN( "6th order exp, alpha & beta vals scaled, and small grid space" ){
      nphon = 6; delta = 0.04; sc = 1.0; scaling = 1.0;
      alpha =  { 0.1, 0.2, 0.4, 0.8, 1.6 };

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      THEN( "contin output matches expected value" ){
        lambda_s = 0.11157823, t_eff = 4.91699518;
        expected = {6.82096404E-5, 7.51470660E-5, 8.20844916E-5, 9.59593429E-5,
         1.23709045E-4, 1.34940719E-4, 1.48666045E-4, 1.62391372E-4, 
          1.89842025E-4, 2.44743331E-4, 2.64063030E-4, 2.90925221E-4, 
          3.17787412E-4, 3.71511794E-4, 4.78960557E-4, 5.05599918E-4, 
          5.57045626E-4, 6.08491333E-4, 7.11382748E-4, 9.17165579E-4, 
          9.26780041E-4, 1.02112863E-3, 1.11547723E-3, 1.30417442E-3, 
          1.68156881E-3};
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff, 1.0e-6 );
      } // THEN
    } // WHEN

    WHEN( "6th order exp, user-defined normalizationand large grid space" ){
      delta = 4.; tbeta = 2.0; sc = 1.0; scaling = 1.0;
      alpha = { 0.1, 0.2, 0.4, 0.8, 1.6 };

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
        alpha, beta, symSab );

      THEN( "contin output matches expected value" ){

        lambda_s = 2.179428E-3; t_eff = 491.1882921;
        expected = {1.37883996E-10, 1.58477481E-10, 1.79070966E-10, 
          2.20257936E-10, 3.02631876E-10, 2.75707903E-10, 3.16885898E-10, 
          3.58063894E-10, 4.40419886E-10, 6.05131869E-10, 5.51175522E-10, 
          6.33495628E-10, 7.15815735E-10, 8.80455949E-10, 1.20973637E-09, 
          1.10139053E-09, 1.26588730E-09, 1.43038407E-09, 1.75937760E-09, 
          2.41736468E-09, 2.19894405E-09, 2.52736456E-09, 2.85578506E-09, 
          3.51262608E-09, 4.82630810E-09};
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff, 1.0e-6 );
      } // THEN
    } // WHEN
  } // GIVEN 


  GIVEN( "simplified water model (test case #9)" ){
    ntempr = 1; nphon = 100; itemp = 0;
    delta = 0.00255; tbeta = 0.444444; tev = 2.55074596e-2; sc = 0.99186670867058835; 
    scaling = 0.99186670867058835; 
    std::vector<double> alpha { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 
      0.100812, 0.151218, 0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 
      0.453654, 0.504060, 0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 
      0.889061, 0.976435, 1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 
      1.707750, 1.873790, 2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 
      3.263080, 3.578320, 3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 
      6.217580, 6.816500, 7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 
      11.83030, 12.96740, 14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 
      22.49220, 24.65260, 27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 
      42.74530, 46.85030, 50.0 },
    beta { 0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, 
      .0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 
      0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 
      1.129090, 1.209740, 1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 
      1.693640, 1.774290, 1.854940, 1.935590, 2.016240, 2.096890, 2.177540, 
      2.258190, 2.338840, 2.419490, 2.500140, 2.580790, 2.669500, 2.767090, 
      2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 3.595360, 3.785490, 
      3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 5.769970, 
      6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 
      10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 17.17330, 
      18.72180, 20.42450, 22.29760, 24.35720, 25.0 },
    rho { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
      0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 0.0686, 
      0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 0.1218, 0.1195, 
      0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 0.0839, 0.0807, 0.07798, 
      0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 0.06652, 0.065, 0.0634, 
      0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 0.05462, 0.0535, 0.0525, 
      0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 0.0459, 0.04478, 0.04366, 
      0.04288, 0.04244, 0.042, 0.0 };
    WHEN( "  " ){

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );
      std::vector<double> correct_alpha_0_to_49_beta_0 { 5.70056684e-4, 
      5.71888708e-4, 5.73720731e-4, 5.77384778e-4, 5.81048825e-4, 5.84712872e-4, 
      5.88951672e-4, 5.93233434e-4, 5.42056424e-4, 4.28244768e-4, 3.06057223e-4, 
      2.95212184e-4, 3.02663494e-4, 2.92398496e-4, 3.05536548e-4, 3.19025735e-4, 
      3.26220873e-4, 3.36551155e-4, 3.49179421e-4, 3.57425983e-4, 3.68816500e-4, 
      3.83479183e-4, 3.94402857e-4, 4.06891766e-4, 4.20261279e-4, 4.32403821e-4, 
      4.42685821e-4, 4.51880855e-4, 4.59388249e-4, 4.67775528e-4, 4.76936118e-4, 
      4.86243516e-4, 4.94620894e-4, 5.09851031e-4, 5.23765935e-4, 5.36531286e-4, 
      5.48284769e-4, 5.59147457e-4, 5.69019309e-4, 5.68792651e-4, 5.58615581e-4, 
      5.39897095e-4, 5.12131779e-4, 4.64702762e-4, 4.11176570e-4, 3.61173172e-4, 
      3.17924984e-4, 2.79777761e-4, 2.46493592e-4, 2.18291694e-4 }, 
      correct_alpha_0_to_49_beta_29 { 7.68099243e-2, 7.70496243e-2, 
      7.72893244e-2, 7.77687244e-2, 7.82481245e-2, 7.87275245e-2, 7.92821246e-2, 
      7.98423458e-2, 7.34097422e-2, 5.90448823e-2, 4.36514932e-2, 4.23617235e-2, 
      4.34464736e-2, 4.23196096e-2, 4.41562638e-2, 4.60477871e-2, 4.71519543e-2, 
      4.86511099e-2, 5.04440988e-2, 5.16881615e-2, 5.33320390e-2, 5.53918648e-2, 
      5.69835987e-2, 5.87731321e-2, 6.06744252e-2, 6.24198691e-2, 6.39279453e-2, 
      6.52948011e-2, 6.64468553e-2, 6.77102160e-2, 6.90717811e-2, 7.04534180e-2, 
      7.17207122e-2, 7.38571845e-2, 7.58283127e-2, 7.76512733e-2, 7.93411960e-2, 
      8.09124408e-2, 8.23410247e-2, 8.24717121e-2, 8.13342100e-2, 7.90998261e-2, 
      7.56977801e-2, 6.97809264e-2, 6.30729464e-2, 5.68179991e-2, 5.14305358e-2, 
      4.66982829e-2, 4.25772650e-2, 3.90777534e-2 }, 
      correct_alpha_49_to_74_beta_64 { 5.01200422e-3, 5.39761518e-3, 
      5.84469397e-3, 6.36259633e-3, 6.96450163e-3, 7.66484981e-3, 8.48120780e-3, 
      9.43557399e-3, 1.05492961e-2, 1.18545543e-2, 1.33788191e-2, 1.51495464e-2, 
      1.71949379e-2, 1.95388039e-2, 2.21987489e-2, 2.51587327e-2, 2.83866371e-2, 
      3.17936243e-2, 3.52301989e-2, 3.84774693e-2, 4.12313570e-2, 4.31193070e-2, 
      4.37377893e-2, 4.27276040e-2, 3.98732704e-2, 3.86829542e-2 };

      THEN( "contin output matches expected value" ){
        lambda_s =   0.23520650571218535; t_eff = 1.9344846581861184;
        REQUIRE( lambda_s == Approx(std::get<0>(output)).epsilon(1e-6) );
        REQUIRE( t_eff    == Approx(std::get<1>(output)).epsilon(1e-6) );
        for ( size_t i = 0; i < correct_alpha_0_to_49_beta_0.size(); ++i ){ 
          REQUIRE( correct_alpha_0_to_49_beta_0[i] == Approx(symSab(0,i,0)).epsilon(1e-6) );
        }
        for ( size_t i = 0; i < correct_alpha_0_to_49_beta_0.size(); ++i ){ 
          REQUIRE( correct_alpha_0_to_49_beta_0[i] == Approx(symSab(0,i,0)).epsilon(1e-6) );
        }
        for ( size_t i = 0; i < correct_alpha_49_to_74_beta_64.size(); ++i ){ 
          REQUIRE( correct_alpha_49_to_74_beta_64[i] == Approx(symSab(64,49+i,0)).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN
  } // GIVEN 

  GIVEN( "extremely simplified water model (test case #9) (simplifiedLeaprInput)" ){
    ntempr = 1; nphon = 6; itemp = 0;
    delta = 0.00255; tbeta = 0.444444; tev = 2.55074596e-2; sc = 0.99186670867058835; 
    scaling = 0.99186670867058835; 
    std::vector<double> alpha { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756 };
    std::vector<double> beta { 0.0, 0.006375, 0.01275, 0.0255, 0.03825, 0.051, 0.06575 }; 
    std::vector<double> rho { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165 };
    WHEN( "  " ){
      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      std::vector<double> expected { 4.2543261427413936e-2, 4.2677975670624702e-2, 
        4.2812689913835460e-2, 4.3082118400256997e-2, 4.3351546886678527e-2, 
        4.3620975373100043e-2, 4.3932667151509258e-2, 6.2405605186179013e-2, 
        6.2601782824467045e-2, 6.2797960462755043e-2, 6.3190315739331079e-2, 
        6.3582671015907102e-2, 6.3975026292483125e-2, 6.4428927494796562e-2, 
        0.10177395635934081, 0.10208912917521698, 0.10240430199109316, 
        0.10303464762284552, 0.10366499325459788, 0.10429533888635020, 
        0.10502456226425980, 0.13029305403062777, 0.13069197038636696, 
        0.13109088674210612, 0.13188871945358452, 0.13268655216506289, 
        0.13348438487654118, 0.13440736781727111, 0.18926073679135993, 
        0.18982576724110425, 0.19039079769084855, 0.19152085859033724, 
        0.19265091948982599, 0.19378098038931460, 0.19508830574362507, 
        0.26409668821687571, 0.26485736219135630, 0.26561803616583685, 
        0.26713938411479798, 0.26866073206375923, 0.27018208001272020, 
        0.27194207077720473 };
      THEN( "contin output matches expected value" ){
        lambda_s = 4.2065831071822970; t_eff = 1.0400408284374374;
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff, 1.0e-6 );
      } // THEN
    } // WHEN
  } // GIVEN 
  GIVEN( "extremely simplified water model (test case #9) (simplifiedLeaprInput)" ){
    ntempr = 1; nphon = 30; itemp = 0;
    delta = 0.00255; tbeta = 0.444444; tev = 2.55074596e-2; sc = 0.99186670867058835; 
    scaling = 0.99186670867058835; 
    std::vector<double> alpha { 0.01008, 0.015, 0.0252 },
      beta { 0.000000, 0.006375, 0.012750, 0.025500 },
      rho { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005 };
    WHEN( "  " ){
      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      std::vector<double> expected { 0.21338252777792166, 0.21401421222255254, 
        0.21464589666718339, 0.21590926555644518, 0.30365460687728840, 
        0.30451809601992325, 0.30538158516255781, 0.30710856344782711, 
        0.46547917133239852, 0.46669713613840202, 0.46791510094440553, 
        0.47035103055641264 };

      THEN( "contin output matches expected value" ){
        lambda_s = 14.500767787022335; t_eff = 1.0123228701111786;
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff, 1.0e-5 );
      } // THEN
    } // WHEN
  } // GIVEN 

} // TEST CASE
