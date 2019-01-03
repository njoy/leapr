#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "leapr.cpp"
#include <unsupported/Eigen/CXX11/Tensor>



void checkSab( const std::vector<double>& correctSab,
  const Eigen::Tensor<double,3>& sab ){

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


void checkPartofSab( int aInit, int aFinal, int bInit, int bFinal, 
  const std::vector<double>& correct, const Eigen::Tensor<double,3>& ssm,
  const double& tol = 1e-6 ){
 
  for ( size_t i = 0; i < correct.size(); ++i ){ 
    int a = aInit == aFinal ? aInit : aInit + i;
    int b = bInit == bFinal ? bInit : bInit + i;
    REQUIRE( correct[i] == Approx(ssm(a,b,0)).epsilon(tol) ); 
  } 

}



TEST_CASE( "leapr" ){
  int nphon, ncold, lat;
  double sps, awr, aws, delta, twt, c, tbeta, dka;
  std::vector<double> alpha, beta, temp, rho, oscE, oscW, kappa,
                      ssmDiag(70), ssmChangeA(70), ssmChangeB(70);

  GIVEN( "simple H in H20 input" ) {
    nphon = 100;
    awr    = 0.99917;  ncold = 0; 
    aws   = 1.1; sps = 3.8883; 
    lat   = 1;            
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406 };       
    beta  = { 0.0, 0.006375, 0.01275, 0.0255, 0.03825, 0.051, 0.06575 };
    temp  = { 296.0 };                                          
    delta = 0.00255;                                           
    rho   = { 0.0000, 0.0005, 0.0010, 0.0020, 0.0035, 0.0050, 0.0075, 0.0100, 
      0.0130, 0.0165, 0.0200, 0.0245, 0.0290, 0.0340, 0.0395, 0.0450, 0.0506, 
      0.0562, 0.0622, 0.0686, 0.0750, 0.0830, 0.0910, 0.0990, 0.1070, 0.1150, 
      0.1197, 0.1214, 0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.0954, 0.0912, 
      0.0871, 0.0839, 0.0807, 0.0779, 0.0757, 0.0735, 0.0716, 0.0697, 0.0680, 
      0.0665, 0.0650, 0.0634, 0.0618, 0.0602, 0.0586, 0.0571, 0.0558, 0.0546, 
      0.0535, 0.0525, 0.0515, 0.0504, 0.0493, 0.0482, 0.0470, 0.0459, 0.0447, 
      0.0436, 0.0428, 0.0424, 0.0420, 0.0000 };
    twt    = 0.055556;     c      = 0.0;    tbeta = 0.444444;   
    oscE   = { 0.205,    0.48};                                 
    oscW   = { 0.166667, 0.333333 };                            
    dka    = 0.01;                            
    kappa  = { 0.1, 0.2, 0.4, 0.7 };                            

    WHEN( "Energy grid is not specified (uniform grid assumed)" ){
      auto out = leapr( nphon, awr, 
          ncold, aws, lat, alpha, beta, 
          temp, delta, rho, twt, c, tbeta, oscE, 
          oscW, dka, kappa );
  
      auto ssm = std::get<2>(out);
      std::vector<double> ssmCorrect { 11.935546503364062, 11.733397895550871, 11.152903960781737, 9.0424432881369619, 6.3611753062620506, 3.8620446900714813, 1.8160057120846576, 9.7708377770648926, 9.6752703537208387, 9.3535876750917559, 8.1349490429353253, 6.4278423334271970, 4.6143335017204183, 2.7828528156853025, 7.5169598382557510, 7.4770846554946662, 7.3345014968107360, 6.7684003276843612, 5.9073227108412087, 4.8555605522066960, 3.6047981998984713, 6.5545915569703386, 6.5267386930920859, 6.4366732105195492, 6.0771160333041898, 5.4673797165151754, 4.7258409922821132, 3.7653042995328416, 5.2779210533938592, 5.2629036386269359, 5.2270468610288656, 5.0338875599034454, 4.7140290044004143, 4.2927131891610735, 3.7130577937826925     };


      checkSab( ssmCorrect, ssm );
    } // WHEN

    WHEN( "Energy grid is specified (uniform grid not assumed)" ){
      std::vector<double> energyGrid(rho.size());
      for (size_t i = 0; i < rho.size(); ++i){
        energyGrid[i] = delta*i;
      }
      auto out = leapr( nphon, awr, 
          ncold, aws, lat, alpha, beta, 
          temp, delta, rho, twt, c, tbeta, oscE, 
          oscW, dka, kappa, energyGrid );
  
      auto ssm = std::get<2>(out);
      std::vector<double> ssmCorrect { 11.935546503364062, 11.733397895550871, 11.152903960781737, 9.0424432881369619, 6.3611753062620506, 3.8620446900714813, 1.8160057120846576, 9.7708377770648926, 9.6752703537208387, 9.3535876750917559, 8.1349490429353253, 6.4278423334271970, 4.6143335017204183, 2.7828528156853025, 7.5169598382557510, 7.4770846554946662, 7.3345014968107360, 6.7684003276843612, 5.9073227108412087, 4.8555605522066960, 3.6047981998984713, 6.5545915569703386, 6.5267386930920859, 6.4366732105195492, 6.0771160333041898, 5.4673797165151754, 4.7258409922821132, 3.7653042995328416, 5.2779210533938592, 5.2629036386269359, 5.2270468610288656, 5.0338875599034454, 4.7140290044004143, 4.2927131891610735, 3.7130577937826925     };

      checkSab( ssmCorrect, ssm );
    } // WHEN


  } // GIVEN 

  GIVEN( "H in H2O input (TEST 09)" ) {
    nphon = 100;                   
    awr    = 0.99917;   ncold = 0;
    aws = 15.85316;   sps = 3.8883;
    lat   = 1;                     
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 0.100812, 0.151218,
      0.201624, 0.25203, 0.302436, 0.352842, 0.403248, 0.453654, 0.50406, 
      0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
      1.07213, 1.17708, 1.29211, 1.41822, 1.55633, 1.70775, 1.87379, 2.05566, 
      2.25506, 2.47352, 2.71295, 2.97546, 3.26308, 3.57832, 3.92390, 4.30266, 
      4.71770, 5.17256, 5.67118, 6.21758, 6.81650, 7.47289, 8.19228, 8.98073, 
      9.84489, 10.7919, 11.8303, 12.9674, 14.2145, 15.5815, 17.0796, 18.7208, 
      20.5203, 22.4922, 24.6526, 27.0216, 29.6175, 32.4625, 35.5816, 38.9991, 
      42.7453, 46.8503, 50.0 }; 
    beta = { 0.00000, 0.006375, 0.01275, 0.02550, 0.03825, 0.05100, 0.06575, 
      0.0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 
      0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.04844, 
      1.12909, 1.20974, 1.29039, 1.37104, 1.45169, 1.53234, 1.61299, 1.69364, 
      1.77429, 1.85494, 1.93559, 2.01624, 2.09689, 2.17754, 2.25819, 2.33884, 
      2.41949, 2.50014, 2.58079, 2.66950, 2.76709, 2.87445, 2.99250, 3.12235, 
      3.26530, 3.42247, 3.59536, 3.78549, 3.99467, 4.22473, 4.47787, 4.75631, 
      5.06258, 5.39939, 5.76997, 6.17766, 6.62607, 7.11924, 7.66181, 8.25862, 
      8.91511, 9.63722, 10.4320, 11.3051, 12.2668, 13.3243, 14.4867, 15.0 }; 
    temp   = { 296.0 };                                             
    delta  = 0.00255;                                              
    rho = { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 
      0.013, 0.0165, 0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 
      0.0622, 0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
      0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 0.0839, 
      0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 0.06652, 0.065, 
      0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 0.05462, 0.0535, 0.0525, 
      0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 0.0459, 0.04478, 0.04366, 
      0.04288, 0.04244, 0.042, 0.0 };          
    twt = 0.055556; c = 0.0; tbeta = 0.444444; 
    oscE = { 0.205, 0.48 };                    
    oscW = { 0.166667, 0.333333 };             
    dka = 0.0;                        
    kappa = { };                            

    auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                      rho, twt, c, tbeta, oscE, oscW, dka, kappa );

    double lambda_s = std::get<0>(out);
    double t_eff = std::get<1>(out);
    auto ssm = std::get<2>(out);

    ssmDiag = { 11.93561, 9.675395, 7.334811, 6.077737, 4.715214, 3.763573, 
    3.135277, 2.543597, 1.945134, 1.529973, 0.961759, 0.585923, 0.352446, 
    0.212400, 0.130001, 8.253134e-2, 5.688670e-2, 4.345888e-2, 3.682223e-2, 
    3.40679e-2, 3.36325e-2, 3.46283e-2, 3.65466e-2, 3.91319e-2, 4.21916e-2, 
    4.56340e-2, 4.93880e-2, 5.33684e-2, 5.75277e-2, 6.17923e-2, 6.60659e-2, 
    7.02334e-2, 7.42183e-2, 7.78721e-2, 8.11006e-2, 8.38128e-2, 8.59146e-2, 
    8.73652e-2, 8.81012e-2, 8.81156e-2, 8.73487e-2, 8.58047e-2, 8.34929e-2, 
    8.04708e-2, 7.67892e-2, 7.25609e-2, 6.78795e-2, 6.28789e-2, 5.76449e-2, 
    5.23256e-2, 4.69968e-2, 4.17509e-2, 3.66860e-2, 3.18639e-2, 2.73362e-2, 
    2.31459e-2, 1.93246e-2, 1.58895e-2, 1.28531e-2, 1.02142e-2, 7.96129e-3, 
    6.07774e-3, 4.53555e-3, 3.30380e-3, 2.73618e-3 };
    ssmChangeB = { 0.507543, 0.508714, 0.509901, 0.512080, 0.514362, 
    0.516854, 0.519580, 0.520541, 0.522258, 0.520150, 0.506676, 0.481421, 
    0.446141, 0.402707, 0.355004, 0.306458, 0.259873, 0.217086, 0.179315, 
    0.147340, 0.121591, 0.101728, 8.70041e-2, 7.65639e-2, 6.95008e-2, 
    6.51034e-2, 6.26512e-2, 6.15678e-2, 6.13886e-2, 6.17923e-2, 6.25461e-2, 
    6.34764e-2, 6.44859e-2, 6.54697e-2, 6.63856e-2, 6.71735e-2, 6.77863e-2, 
    6.82127e-2, 6.83930e-2, 6.83276e-2 };
    ssmChangeA = { 5.5741e-4, 8.28335e-4, 1.38765e-3, 1.81319e-3, 2.75597e-3,
    4.10379e-3, 5.43260e-3, 8.03115e-3, 1.05535e-2, 1.30041e-2, 1.53847e-2, 
    1.77000e-2, 1.99510e-2, 2.21390e-2, 2.42686e-2, 2.63383e-2, 2.85428e-2, 
    3.08854e-2, 3.33695e-2, 3.59893e-2, 3.87456e-2, 4.16348e-2, 4.46483e-2, 
    4.77815e-2, 5.10161e-2, 5.43383e-2, 5.77140e-2, 6.11324e-2, 6.45536e-2, 
    6.79373e-2, 7.12310e-2, 7.43908e-2, 7.73717e-2, 8.01064e-2, 8.25233e-2, 
    8.45724e-2, 8.61959e-2, 8.73502e-2, 8.79674e-2, 8.79771e-2 };
     
    checkPartofSab(  0, 65,  0, 65, ssmDiag,    ssm, 5e-5 );
    checkPartofSab( 29, 29,  0, 40, ssmChangeB, ssm, 1e-6 );
    checkPartofSab(  0, 40, 40, 40, ssmChangeA, ssm, 1e-6 );
   
    REQUIRE( 0.23520650 == Approx(lambda_s).epsilon(1e-6) );        
    REQUIRE( 1.93448465 == Approx(t_eff).epsilon(1e-6) );        

    
    sps += 1.0; 
  } // GIVEN 

  GIVEN( "Para hydrogen at 20K (TEST 22)" ) {
    nphon = 100;                 
    double smin = 2.0e-38;                                          
    awr = 0.99917; ncold = 2;
    aws = 0.0; sps = 0.0;         
    lat = 0;                      
    alpha = { 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 
    0.7, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 
    60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 
    360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600 };
    beta = { 0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 2.5, 
    3, 3.5, 4, 4.5, 5, 5.5, 5.7, 5.8, 5.9, 6, 6.25, 6.5, 7, 7.5, 8, 8.4, 8.5, 
    8.6, 8.7027, 8.8, 8.9, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 
    14.5, 15, 15.5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 25.7, 26.1, 26.5, 
    27, 28, 29, 30, 35, 37.5, 40, 45.5, 50, 52.2, 55, 60, 65, 70, 75, 80, 85, 
    90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 160, 170, 
    180, 190, 200, 210, 220, 230, 240, 250, 260, 280, 300 };
    temp   = { 20.0 };                                             
    delta  = 0.00025;                                             
    rho = { 0.0, 0.01563, 0.0625, 0.141, 0.25, 0.391, 0.5625, 0.766, 1., 1.266, 
    1.5625, 1.89, 2.25, 2.64, 3.0625, 3.52, 4., 4.6, 5.5, 7.0, 8.5, 9.2, 9.5, 
    9.4, 9.2, 8.9, 8.5, 8.0, 7.5, 7.05, 6.7, 6.4, 6.2, 6.1, 6.2, 6.45, 6.7, 
    6.95, 7.1, 6.55, 5.5, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };     
    twt = 0.025; c = 40.0; tbeta = 0.475;      
    oscE = { };                                
    oscW = { };                                
    dka = 0.05;                     
    kappa = { 9.29815e-2, 9.29815e-2, 9.31478e-2, 9.43837e-2, 9.5458e-2, 
    9.68871e-2, 9.86644e-2, 0.100813, 0.103358, 0.106334, 0.109777, 0.11373, 
    0.118249, 0.123395, 0.129241, 0.135875, 0.143399, 0.151932, 0.161618, 
    0.172622, 0.185143, 0.199416, 0.215717, 0.234376, 0.255779, 0.280385, 
    0.308731, 0.341443, 0.37924, 0.422934, 0.473417, 0.53162, 0.598441, 0.674614, 
    0.79867, 1.118, 1.4087, 1.6633, 1.8799, 2.0554, 2.1889, 2.2803, 2.3306, 
    2.3418, 2.3164, 2.2581, 2.1706, 2.0585, 1.9264, 1.779, 1.6213, 1.4581, 
    1.2938, 1.1327, 0.97879, 0.83545, 0.70565, 0.5918, 0.49576, 0.41882, 0.36172, 
    0.32464, 0.30727, 0.30884, 0.32812, 0.36355, 0.41325, 0.4751, 0.54681, 
    0.62601, 0.71023, 0.79708, 0.88421, 0.96942, 1.0507, 1.1261, 1.1943, 1.2538, 
    1.3037, 1.3432, 1.3719, 1.3897, 1.3968, 1.3935, 1.3805, 1.3586, 1.3289, 
    1.2924, 1.2505, 1.2045, 1.1557, 1.1055, 1.0551, 1.0059, 0.95899, 0.91541, 
    0.87608, 0.84174, 0.81298, 0.79025, 0.7738, 0.76372, 0.75995, 0.76226, 
    0.77032, 0.78364, 0.80164, 0.82367, 0.84898, 0.87681, 0.90637, 0.93684, 
    0.96745, 0.99744, 1.0261, 1.0529, 1.0771, 1.0984, 1.1163, 1.1306, 1.141, 
    1.1476, 1.1504, 1.1493, 1.1447, 1.1369, 1.126, 1.1126, 1.0971, 1.0799, 
    1.0616, 1.0426, 1.0233, 1.0044, 0.98614, 0.96901, 0.95336, 0.93949, 0.92766, 
    0.91806, 0.91082, 0.90602, 0.90365, 0.90367, 0.90598, 0.91043, 0.91682, 
    0.92491, 0.93444, 0.94513, 0.95667, 0.96876, 0.98109, 0.99335, 1.0053, 
    1.0166, 1.027, 1.0364, 1.0445, 1.0513, 1.0566, 1.0603, 1.0625, 1.0631, 
    1.0622, 1.0599, 1.0563, 1.0515, 1.0457, 1.0391, 1.0318, 1.0241, 1.0161, 
    1.008, 1.0001, 0.9925, 0.98538, 0.97887, 0.97311, 0.9682, 0.96422, 0.96122, 
    0.95924, 0.95828, 0.95832, 0.95932, 0.96122, 0.96394, 0.96739, 0.97145, 
    0.97602, 0.98096, 0.98614, 0.99145, 0.99674, 1.0019, 1.0068, 1.0114, 1.0155, 
    1.0192, 1.0222, 1.0246, 1.0264, 1.0275, 1.0279, 1.0277, 1.0268, 1.0254, 
    1.0234, 1.021, 1.0182, 1.0151, 1.0117, 1.0082, 1.0047, 1.0012, 0.99774, 
    0.9945, 0.99151, 0.98882, 0.98648, 0.98454, 0.98302, 0.98193, 0.9813, 
    0.98111, 0.98135, 0.9820, 0.98303, 0.98441, 0.98608, 0.98801, 0.99013, 
    0.9924, 0.99475, 0.99713, 0.99949, 1.0018, 1.0039, 1.0059, 1.0077, 1.0092, 
    1.0105, 1.0115, 1.0122, 1.0126, 1.0127, 1.0125, 1.012, 1.0113, 1.0104, 
    1.0092, 1.0079, 1.0065, 1.0049, 1.0034, 1.0017, 1.0001, 0.9986, 0.99715, 
    0.9958, 0.9946, .99355, 0.99269, 0.99201, 0.99154, 0.99126, 0.99119, 0.99132, 
    0.99163, 0.99211, 0.99274, 0.99351, 0.9944, 0.99537, 0.99641, .99749, 
    0.99859, 0.99967, 1.0007, 1.0017, 1.0026, 1.0035, 1.0042, 1.0048, 1.0053, 
    1.0056, 1.0058, 1.0058, 1.0058, 1.0056, 1.0053, 1.0048, 1.0043, 1.0037, 
    1.0031, 1.0024, 1.0016, 1.0009, 1.0001, 0.9994 };     

    auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                      rho, twt, c, tbeta, oscE, oscW, dka, kappa );
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);
    REQUIRE( 0.18256462793803138 == Approx(lambda_s).epsilon(1e-6) );        

    ssmDiag = { 0.6180057, 4.031768e-2, 1.836145e-3, 2.415306e-4, 6.909736e-5, 
    7.438264e-5, 2.214376e-5, 1.185262e-5, 1.269748e-5, 2.048916e-5, 4.175543e-5, 
    7.698098e-5, 1.493513e-4, 3.505935e-4, 7.392398e-4, 1.717140e-3, 4.435776e-3, 
    6.291108e-3, 6.696370e-3, 6.456541e-3, 5.745364e-3, 4.388677e-3, 4.206325e-3, 
    4.726295e-3, 5.695450e-3, 6.316777e-3, 7.815625e-3, 6.644041e-3, 5.012789e-3, 
    1.957844e-3, 6.212740e-4 };
    ssmChangeA = { 3.347467e-2, 1.793655e-2, 3.330896e-3, 4.149734e-4, 
    3.785414e-5, 2.563831e-6, 1.557010e-7, 2.428126e-8, 6.911702e-9, 1.695093e-9, 
    3.424796e-10, 5.975504e-11, 9.382899e-12, 1.358338e-12, 1.839243e-13, 
    2.351306e-14, 2.855355e-15, 3.307825e-16, 3.672717e-17, 3.936153e-18, 
    4.132063e-19, 4.376742e-20, 4.915711e-21, 6.184114e-22, 8.881905e-23, 
    1.406363e-23, 2.324160e-24, 3.849762e-25, 6.260000e-26 };
    ssmChangeB = { 1.070048e-2, 7.926089e-3, 5.508271e-3, 3.649397e-3, 
    4.022429e-3, 5.014804e-3, 5.265416e-3, 4.447225e-3, 3.000786e-3, 1.777324e-3, 
    8.733902e-4, 3.901252e-4, 1.601165e-4, 6.258589e-5, 2.543928e-5, 1.084678e-5, 
    4.893507e-6, 3.277546e-6, 7.217017e-6, 1.737234e-5, 2.911928e-5, 3.524165e-5, 
    2.368799e-5, 7.715643e-6, 1.448075e-6, 1.784518e-7, 4.546030e-8, 2.525618e-8, 
    8.649551e-9, 1.695710e-9, 1.397482e-9, 2.053214e-8, 4.782204e-8, 1.210798e-8, 
    2.801668e-10 };

    checkPartofSab(  0, 31,  0, 31,  ssmDiag,    ssm, 1e-6 );
    checkPartofSab( 30, 59, 58, 58,  ssmChangeA, ssm, 1e-6 );
    checkPartofSab( 30, 30, 70, 105, ssmChangeB, ssm, 1e-6 );
  } // GIVEN 

  GIVEN( "Beryllium Metal (pg 25 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
   "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
   "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {
    nphon = 150;                 
    awr    = 8.93478; ncold = 0;
    aws = 0.0; sps = 0.0;         
    lat = 1;                      
    alpha = { 0.01008, 0.01500, 0.02520, 0.03300, 0.05040, 0.07560, 0.10080, 
              0.15000, 0.25203, 0.33000, 0.50406, 0.75609, 1.00812, 1.26015, 
              1.51218, 1.76421, 2.01624, 2.26827, 2.52030, 2.77233, 3.02436, 
              3.28163, 3.54445, 3.81261, 4.08651, 4.36617, 4.65167, 4.94311, 
              5.24081, 5.54466, 5.85486, 6.17161, 6.49501, 6.82527, 7.16239, 
              7.50666, 7.85809, 8.21688, 8.58323, 8.95735, 9.33922, 9.72916, 
              10.1276, 10.5338, 10.9492, 11.3726, 11.8051, 12.2466, 12.9830, 
              13.1580, 13.6288, 14.1086, 14.5986, 15.0986, 15.6097, 16.1309, 
              16.6642, 17.2076, 17.7631, 18.3, 19.0, 20.0, 21.0, 22.0, 23.0, 
              24.0, 25.0, 26.0, 27.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 
              42.0, 44.0, 46.0, 48.0, 50.0, 52.5, 55.0, 60.0, 62.5, 65.0, 67.5, 
              70.0, 75.0, 80.0 };                               
    beta = { 0.00000, 1.51218e-1, 3.02436e-1, 4.53654e-1, 6.04872e-1, 
            7.56090e-1, 9.07308e-1, 1.05853, 1.20974, 1.36096, 1.51218, 1.66340, 
            1.81462, 1.890225, 1.96583, 2.04144, 2.11705, 2.19266, 2.26827, 
            2.34388, 2.41949, 2.49510, 2.57071, 2.646315, 2.72192, 2.79753, 
            2.87314, 2.94875, 3.02436, 3.17558, 3.32679, 3.47801, 3.62923, 
            3.78045, 3.93167, 4.08288, 4.24146, 4.40790, 4.58251, 4.76558, 
            4.95763, 5.15915, 5.37046, 5.48135, 5.59224, 5.708525, 5.82481, 
            6.06878, 6.32474, 6.59320, 6.87487, 7.17025, 7.48015, 7.80527, 
            8.14631, 8.50399, 8.87932, 9.27289, 9.68581, 1.01195e1, 1.05732e1, 
            1.10500e1, 1.15500e1, 1.20742e1, 1.26247e1, 1.32023e1, 1.38072e1, 
            1.44423e1, 1.51087e1, 1.58073e1, 1.65412e1, 1.73094e1, 1.81169e1, 
            1.89627e1, 1.98509e1, 2.07824e1, 2.17592e1, 2.27835e1, 2.38581e1, 
            2.49862e1, 2.61688e1, 2.74098e1, 2.87112e1, 3.00772e1, 3.15088e1, 
            3.30119e1, 3.45876e1, 3.62409e1, 3.79749e1, 3.97945e1, 42.0, 44.0, 
            46.0, 48.0, 50.0, 52.5, 55.0, 57.5, 60.0, 65.0, 70.0, 75.0, 80.0 };
    temp   = { 296.0 };                                             
    delta  = 0.0015361;                                            
    rho = { 0.0, 0.0062, 0.025, 0.06, 0.115, 0.17, 0.235, 0.308, 0.39, 0.475, 
            0.595, 0.775, 1.052, 1.1482, 1.6233, 1.5653, 1.7423, 2.4289, 2.8072, 
            3.2863, 3.7577, 4.4397, 5.4924, 6.3315, 7.6421, 9.5339, 12.1016, 
            15.0553, 22.0154, 26.6382, 29.3387, 32.6036, 35.777, 37.6536, 
            40.3845, 37.1196, 32.39, 29.8269, 21.2831, 18.9488, 13.0186, 
            7.33848, 12.5075, 18.9031, 23.5564, 25.6923, 25.5855, 43.0544, 
            17.8656, 6.81975, 3.98963, 2.99184, 1.6843, 0.65772, 0.0 }; 
    twt = 0.0; c = 0.0; tbeta = 1.0;      
                                    
    oscE = { };                                
    oscW = { };                                
    dka = 0.0;                     
    kappa = { };     

    auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                      rho, twt, c, tbeta, oscE, oscW, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    ssmDiag = { 4.359730e-4, 7.268573e-4, 1.378634e-3, 1.758590e-3, 2.746460e-3, 
    5.160857e-3, 7.213335e-3, 1.341711e-2, 2.482782e-2, 4.035686e-2, 7.538107e-2, 
    0.1562404, 0.2523771, 0.2874854, 0.3071275, 0.3131106, 0.2826244, 0.2321774, 
    0.1898775, 0.1473543, 0.1172308 };
    ssmChangeB = { 5.251788e-9, 5.758381e-9, 6.341055e-9, 6.668745e-9, 
    7.012559e-9, 7.391516e-9, 7.790471e-9, 8.695648e-9, 9.754232e-9, 1.099733e-8, 
    1.246419e-8, 1.420376e-8, 1.627976e-8, 1.876997e-8, 2.177230e-8, 2.541455e-8, 
    2.986444e-8, 3.532792e-8, 4.208514e-8, 5.050687e-8, 6.103243e-8, 7.434134e-8, 
    9.126652e-8, 1.129272e-7, 1.409181e-7, 1.773543e-7, 2.250379e-7, 2.881027e-7, 
    3.721754e-7, 4.850630e-7 };
    ssmChangeA = { 5.762397e-22, 2.405963e-21, 9.553490e-21, 3.621451e-20, 
    1.313920e-19, 4.578065e-19, 1.534946e-18, 4.964381e-18, 1.552183e-17, 
    4.700931e-17, 1.380635e-16, 3.939845e-16, 1.094604e-15, 2.957767e-15, 
    7.806162e-15, 2.007872e-14, 5.048931e-14, 1.241469e-13, 5.110596e-13, 
    7.049273e-13, 1.631523e-12, 3.697868e-12, 8.227070e-12, 1.796486e-11, 
    3.857360e-11, 8.133015e-11, 1.688790e-10, 3.444321e-10, 6.918538e-10, 
    1.320041e-9, 2.948102e-9, 8.659421e-9, 2.358561e-8, 5.999104e-8, 1.433741e-7, 
    3.236886e-7, 6.935889e-7, 1.416451e-6, 2.767136e-6, 5.188240e-6, 1.631094e-5 };

    checkPartofSab(  0, 21,   0,  21, ssmDiag,    ssm, 1e-6 );
    checkPartofSab( 30, 71, 100, 100, ssmChangeA, ssm, 1e-6 );
    checkPartofSab( 89, 89,  40,  70, ssmChangeB, ssm, 1e-6 );
 
    REQUIRE( 0.73178450430084030 == Approx(lambda_s).epsilon(1e-6) );        
  } // GIVEN 

  GIVEN( "Beryllium Oxide (pg 39 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
    "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
    "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {
    nphon = 100;                 
    awr    = 8.93478; ncold = 0;
    aws = 15.858; sps = 3.7481; int mss = 1;         
    lat = 1;                      
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.0504, 0.0756, 0.1008, 0.15, 0.252, 
    0.33, 0.504, 0.756, 1.008, 1.260, 1.512, 1.764, 2.016, 2.268, 2.520, 2.772, 
    3.024, 3.282, 3.544, 3.813, 4.087, 4.366, 4.652, 4.943, 5.241, 5.545, 5.855, 
    6.172, 6.495, 6.825, 7.162, 7.507, 7.858, 8.217, 8.583, 8.957, 9.339, 9.729, 
    10.13, 10.53, 10.95, 11.37, 11.81, 12.25, 12.69, 13.16, 13.63, 14.11, 14.60, 
    15.10, 15.61, 16.13, 16.66, 17.21, 17.76, 18.3, 19, 20, 21, 22, 23, 24, 25, 
    26, 27, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52.5, 55, 57.5, 60, 
    62.5, 65, 70, 75, 80 };                               
    beta = { 0, 0.1513, 0.3025, 0.4537, 0.6049, 0.7561, 0.9073, 1.059, 1.210, 
    1.361, 1.512, 1.5875, 1.663, 1.7390, 1.815, 1.8905, 1.966, 2.0415, 2.117, 
    2.1925, 2.268, 2.3435, 2.419, 2.4950, 2.571, 2.6465, 2.722, 2.7975, 2.873, 
    2.9485, 3.024, 3.1000, 3.176, 3.2515, 3.327, 3.4025, 3.478, 3.5850, 3.629, 
    3.7045, 3.780, 3.8560, 3.932, 4.0075, 4.083, 4.1620, 4.241, 4.3245, 4.408, 
    4.4955, 4.583, 4.6745, 4.766, 4.8620, 4.958, 5.0585, 5.159, 5.371, 5.592, 
    5.825, 6.069, 6.325, 6.593, 6.875, 7.170, 7.480, 7.805, 8.146, 8.504, 8.879, 
    9.273, 9.686, 10.12, 10.57, 11.05, 11.55, 12.07, 12.62, 13.20, 13.81, 14.44, 
    15.11, 15.81, 16.54, 17.31, 18.12, 18.96, 19.85, 20.78, 21.76, 22.78, 23.86, 
    24.99, 26.17, 27.41, 28.71, 30.08, 31.51, 33.01, 34.59, 36.24, 37.98, 39.80, 
    42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 65, 70, 75, 80 };
    temp   = { 296.0 };                                             
    delta  = 0.0016518;                                            
    rho = { 0.0, 0.3, 0.7, 0.9, 1.0, 1.2, 1.6, 2.0, 2.2, 3.0, 3.5, 4.5, 5.5, 
    6.8, 8.0, 9.2, 10.9, 12.9, 15.5, 18.6, 22.0, 26.0, 30.5, 35.0, 39.0, 40.0, 
    34.0, 28.0, 26.0, 24.4, 23.0, 21.3, 19.8, 17.0, 14.1, 12.0, 10.0, 9.0, 9.0, 
    8.5, 7.5, 6.0, 4.6, 3.1, 1.6, 0.5, 0.0, 0.0, 4.0, 15.0, 38.0, 52.0, 70.0, 
    105.0, 165.0, 230.0, 200.0, 170.0, 145.0, 136.0, 134.0, 112.0, 96.0, 89.0, 
    84.0, 75.0, 87.0, 81.0, 66.0, 59.0, 68.0, 105.0, 95.0, 97.0, 135.0, 163.0, 
    130.0, 111.0, 92.0, 67.0, 45.0, 19.0, 7.0, 0.0 }; 
    twt = 0.0; c = 0.0; tbeta = 1.0;      
                                       
    oscE = { };                                
    oscW = { };                                
    dka = 0.0;                     
    kappa = { };     

    auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                      rho, twt, c, tbeta, oscE, oscW, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    ssmDiag = { 2.79995e-3, 2.24993e-3, 1.43145e-3, 1.54553e-3, 2.23734e-3, 
    3.68195e-3, 5.59778e-3, 9.20979e-3, 1.83318e-2, 2.89240e-2, 4.97225e-2, 
    6.98011e-2, 7.59658e-2, 6.84755e-2, 6.40908e-2, 6.21202e-2, 6.01752e-2, 
    5.75049e-2, 5.41214e-2, 4.93725e-2, 4.58363e-2 };
    ssmChangeB = { 4.28695e-5, 5.70149e-5, 7.61004e-5, 1.01837e-4, 1.36809e-4, 
    1.84266e-4, 2.49036e-4, 3.36505e-4, 4.54908e-4, 6.15442e-4, 8.30266e-4, 
    1.11891e-3, 1.50079e-3, 2.08894e-3, 2.75954e-3, 3.57197e-3, 4.53361e-3, 
    5.64596e-3, 6.90341e-3, 8.29236e-3, 9.79100e-3, 1.13694e-2, 1.29906e-2, 
    1.69426e-2, 2.01680e-2, 2.20428e-2, 2.22372e-2 };
    ssmChangeA = { 7.83271e-22, 3.42838e-21, 1.36709e-20, 5.09538e-20, 
    1.77115e-19, 5.78380e-19, 1.79913e-18, 5.31309e-18, 1.50672e-17, 4.10182e-17, 
    1.07524e-16, 2.72927e-16, 6.70289e-16, 1.60035e-15, 3.72035e-15, 8.45264e-15, 
    1.87041e-14, 4.05420e-14, 8.59741e-14, 1.78894e-13, 3.6550e-13 };

    checkPartofSab(  0, 21,   0,  21, ssmDiag,    ssm, 1e-6 );
    checkPartofSab( 20, 41, 116, 116, ssmChangeA, ssm, 1e-6 );
    checkPartofSab( 89, 89,  90, 117, ssmChangeB, ssm, 1e-6 );
 
    REQUIRE( 0.49335785747954103 == Approx(lambda_s).epsilon(1e-6) );        
  } // GIVEN 

  GIVEN( "Graphite (pg 51 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
    "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
    "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {
    nphon = 100;                 
    awr    = 11.898; ncold = 0;
    aws = 0; sps = 0; int mss = 0;         
    lat = 1;                      
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.0504, 0.0756, 0.1008, 0.15, 
    0.25203, 0.33, 0.50406, 0.75609, 1.00812, 1.26015, 1.51218, 1.76421, 
    2.01624, 2.27331, 2.53552, 2.80297, 3.07577, 3.35401, 3.6379, 3.92733, 
    4.22271, 4.52383, 4.83111, 5.14443, 5.46411, 5.79013, 6.12261, 6.46185, 
    6.80783, 7.16077, 7.52067, 7.88783, 8.26234, 8.64432, 9.03396, 9.43136, 
    9.83673, 10.2506, 10.6719, 11.1024, 11.5409, 11.9886, 12.4452, 12.911, 
    13.3858, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32.5, 35, 37.5, 
    40, 42.5, 45, 47.5, 50, 52.5, 55, 60 };                               
    beta = { 0.000000, 0.100812, 0.201624, 0.302436, 0.403248, 0.504060, 
    0.604872, 0.705684, 0.806496, 0.907307, 1.00812, 1.10893, 1.20974, 1.31055, 
    1.41137, 1.51218, 1.61299, 1.71380, 1.81461, 1.91543, 2.01624, 2.11705, 
    2.21786, 2.31867, 2.41949, 2.52030, 2.62111, 2.72192, 2.82273, 2.92354, 
    3.02436, 3.12517, 3.22598, 3.32679, 3.42760, 3.52842, 3.62923, 3.73004, 
    3.83085, 3.93167, 4.03248, 4.13329, 4.24378, 4.36485, 4.49762, 4.64309, 
    4.80248, 4.97719, 5.16873, 5.37862, 5.60867, 5.73473, 5.86080, 5.99896, 
    6.13713, 6.28855, 6.43997, 6.60591, 6.77184, 6.95376, 7.13567, 7.33502, 
    7.53438, 7.75289, 7.97140, 8.21088, 8.45036, 8.97529, 9.55052, 10.1810, 
    10.8726, 11.6297, 12.4593, 13.3697, 14.3667, 15.4595, 16.6571, 17.9697, 
    19.4093, 20.9860, 22.7139, 24.6082, 26.6849, 28.9602, 31.4533, 34.1873, 
    37.1825, 40.4659, 45, 50, 55, 60, 65, 70, 75, 80 };
    temp   = { 296.0 };                                 
    delta  = 0.005485;                                 
    rho = { 0, 0.346613, 1.4135, 3.03321, 3.25901, 3.38468, 3.48269, 3.76397, 
    4.05025, 4.84696, 7.35744, 5.88224, 4.63255, 4.48287, 5.80642, 4.63802, 
    4.28503, 3.92079, 4.91352, 5.53836, 7.51076, 5.31651, 5.40525, 5.20376, 
    5.3276, 7.17251, 3.31813, 4.50126, 5.04663, 4.2089, 2.91985, 4.65109, 
    13.1324, 7.25016, 6.5662, 5.47181, 5.06137, 5.19813, 0.457086, 0 }; 
    twt = 0.0; c = 0.0; tbeta = 1.0;   
                               
    oscE = { };                        
    oscW = { };                        
    dka = 0.0;                
    kappa = { };                       

    auto out = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                      rho, twt, c, tbeta, oscE, oscW, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    ssmDiag = { 1.897534e-3, 2.962635e-3, 5.195673e-3, 7.171099e-3, 1.147201e-2, 
    1.755189e-2, 2.364685e-2, 3.176487e-2, 4.241821e-2, 4.512548e-2, 5.639733e-2, 
    6.788740e-2, 7.598406e-2, 8.040455e-2, 8.505987e-2, 8.774893e-2, 8.989972e-2, 
    9.100598e-2, 9.311458e-2, 9.473819e-2 };
    ssmChangeB = { 2.506645e-5, 2.696827e-5, 2.899493e-5, 3.136419e-5, 
    3.390166e-5, 3.688178e-5, 4.008670e-5, 4.802162e-5, 5.819702e-5, 7.147976e-5, 
    8.900321e-5, 1.122962e-4, 1.435889e-4, 1.861373e-4, 2.444525e-4, 3.251682e-4, 
    4.376761e-4, 5.956025e-4, 8.185592e-4, 1.133181e-3 };
    ssmChangeA = { 5.072964e-13, 1.001498e-12, 1.932624e-12, 3.652261e-12, 
    6.767699e-12, 1.231643e-11, 2.203782e-11, 3.880941e-11, 6.733347e-11, 
    1.151841e-10, 1.944429e-10, 3.242862e-10, 5.338676e-10, 8.697364e-10, 
    1.4006873e-9, 2.2335475e-9, 3.5262728e-9, 5.5148241e-9, 8.5437866e-9, 
    1.4678455e-8 };

    checkPartofSab(  0, 21,  0, 21, ssmDiag,    ssm, 1e-6 );
    checkPartofSab( 30, 51, 95, 95, ssmChangeA, ssm, 1e-6 );
    checkPartofSab( 71, 71, 60, 81, ssmChangeB, ssm, 1e-6 );
 
    REQUIRE( 0.6675520 == Approx(lambda_s).epsilon(1e-6) ); 
  } // GIVEN 
  
  GIVEN( "H in ZrH (pg 75 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
    "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
    "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {
    nphon = 100;                 
    awr    = 0.99917; ncold = 0;   
    aws = 0; sps = 0; int mss = 0; 
    lat = 1;                       
    alpha = { 0.504060, 1.00812, 1.51218, 2.01624, 2.52030, 3.02436, 3.52842, 
    4.03248, 4.53654, 5.04060, 5.60867, 6.24893, 6.97044, 7.78359, 8.69997, 
    9.73279, 10.8968, 12.2083, 13.6872, 15.3526, 17.2308, 19.3468, 21.7320, 
    24.4197, 27.4491, 30.8636, 34.7106, 39.0465, 43.9338, 49.4412, 55.6482, 
    62.6425, 70.5260, 79.4105, 89.4242, 100.708, 113.423, 127.759, 143.909, 
    162.116, 170, 180, 190, 200, 210, 220, 230, 240 };                               
    beta = { 0.0, 0.0790678, 0.158134, 0.237200, 0.316277, 0.395344, 0.474411, 
    0.513939, 0.553478, 0.593016, 0.632545, 0.672083, 0.711611, 0.751150, 
    0.790678, 0.830217, 0.869755, 0.948822, 1.02788, 1.10691, 1.18605, 1.26509, 
    1.34412, 1.46278, 1.58134, 1.69999, 1.81855, 1.93720, 2.05576, 2.17441, 
    2.29297, 2.41162, 2.53018, 2.88594, 3.20229, 3.51854, 3.87430, 4.15103, 
    4.46738, 4.58594, 4.66507, 4.74411, 4.82314, 4.90218, 4.98132, 5.06035, 
    5.13939, 5.21853, 5.29757, 5.33708, 5.37660, 5.41612, 5.45574, 5.53478, 
    5.61381, 5.69295, 5.77199, 5.81150, 5.85102, 5.89054, 5.93016, 5.96968, 
    6.00920, 6.08823, 6.16727, 6.24641, 6.32545, 6.91842, 7.51150, 8.10447, 
    9.09283, 9.40908, 9.72543, 10.0417, 10.3584, 10.6740, 10.7536, 10.8322, 
    10.9119, 11.0691, 11.2274, 11.3857, 11.5440, 11.6226, 11.7023, 11.7809, 
    11.8605, 11.9392, 12.0188, 12.3343, 12.6509, 13.0461, 13.6388, 14.2326, 
    14.8254, 15.4182, 15.7740, 16.0110, 16.2489, 16.4858, 16.7227, 16.9606, 
    17.1975, 17.4344, 17.6713, 17.9092, 18.1855, 18.5020, 18.9768, 19.7672, 
    20.5576, 20.9527, 21.3479, 21.7441, 22.1393, 22.5345, 22.9297, 23.3248, 
    23.7200, 24.5114, 25.3018, 25.6970, 26.0921, 26.4883, 26.8835, 27.2787, 
    27.6739, 28.0691, 28.4642, 28.8594, 29.2556, 29.6508, 30.2436, 30.8363, 
    31.2315, 31.6277, 32.0229, 32.4181, 32.8133, 33.2085, 33.6036, 33.9988, 
    34.3950, 34.7902, 35.5806, 36.3709, 37.1623, 37.9527, 38.7430, 39.5344, 
    40.25, 41, 41.75, 42.5, 43.25, 44, 44.75, 45.5, 46.25, 47, 47.75, 48.5, 
    49.25, 50, 51.75, 52.5, 53.25, 54, 54.75, 55.5, 56.25, 57, 57.75, 58.5, 
    59.25, 60, 60.75, 61.5, 62.25, 63, 63.75, 64.5, 65.25, 66, 66.75, 67.5, 
    68.25, 69, 69.75, 70.5, 71.25, 72, 72.75, 73.5, 74.25, 75, 75.75, 76.5, 
    77.25, 78.0 };
    temp   = { 296.0 };                  
    delta  = 0.001;                     
    rho = { 0, 0.000875, 0.0035, 0.008, 0.015, 0.0235, 0.0340, 0.046, 0.061, 
    0.078, 0.094, 0.116, 0.144, 0.1606, 0.1969, 0.2606, 0.3479, 0.3559, 0.3500, 
    0.3322, 0.3328, 0.2911, 0.1617, 0.1431, 0.1248, 0.09738, 0.06067, 0.1221, 
    0.1495, 0.07219, 0.01443, 0.0001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0499, 2.010, 3.560, 
    4.790, 5.995, 7.250, 8.550, 9.640, 11.91, 13.52, 16.04, 19.79, 26.10, 29.39, 
    30.82, 32.21, 31.75, 33.14, 35.65, 33.34, 36.27, 38.18, 38.75, 39.48, 28.99, 
    23.29, 25.18, 26.59, 27.86, 27.89, 29.44, 25.86, 23.33, 24.66, 27.51, 37.94, 
    60.77, 26.66, 18.54, 14.51, 11.48, 9.53, 7.53, 5.449, 3.838, 8.497, 0 }; 
    twt = 0.0; c = 0.0; tbeta = 1.0;   
                                
    oscE = { };                        
    oscW = { };                        
    dka = 0.0;                
    kappa = { };                       

    auto out1 = leapr( nphon, awr, ncold, aws, lat, alpha, beta, temp, delta, 
                       rho, twt, c, tbeta, oscE, oscW, dka, kappa );
  
    double lambda_s = std::get<0>(out1);
    double t_eff    = std::get<1>(out1);
    auto ssm = std::get<2>(out1);

    ssmDiag = { 6.543001e-3, 1.234760e-2, 1.865458e-2, 2.354496e-2, 2.786202e-2, 
    3.090729e-2, 3.575689e-2, 3.581211e-2, 3.899987e-2, 4.545964e-2, 5.317485e-2, 
    4.810721e-2, 4.163679e-2, 3.456424e-2, 2.974848e-2, 2.240912e-2, 1.181849e-2, 
    7.801597e-3, 4.255764e-3, 5.039730e-3 };
    ssmChangeB = { 2.720975e-3, 4.249577e-3, 5.671485e-3, 6.486683e-3, 
    6.394040e-3, 5.466360e-3, 4.093040e-3, 2.713601e-3, 1.643436e-3, 8.276244e-4, 
    2.209183e-3, 4.008170e-3, 6.482835e-3, 9.189646e-3, 1.138098e-2, 1.236226e-2, 
    1.182379e-2, 1.000693e-2, 7.541887e-3, 5.125047e-3 };
    ssmChangeA = { 6.220088e-7, 2.156774e-6, 7.168986e-6, 2.267933e-5, 
    6.779655e-5, 1.899118e-4, 4.937059e-4, 1.179309e-3, 2.558631e-3, 4.976e-3, 
    8.550100e-3, 1.276913e-2, 1.627568e-2, 1.734380e-2, 1.509877e-2, 1.046399e-2, 
    5.608246e-3, 2.248871e-3, 6.509956e-4, 1.304707e-4 };

    checkPartofSab(  0, 21,  0, 21,   ssmDiag,    ssm, 1e-6 );
    checkPartofSab( 20, 51, 199, 199, ssmChangeA, ssm, 1e-6 );
    checkPartofSab( 30, 30, 110, 131, ssmChangeB, ssm, 1e-6 );
 
    REQUIRE( 0.21630210 == Approx(lambda_s).epsilon(1e-6) ); 
  } // GIVEN 


  GIVEN( "artificial input to test rho energy grid" ) {
    nphon = 100;             
    awr    = 0.99917; ncold = 0;
    aws   = 1.1; sps = 3.8883; 
    lat   = 1;               
    alpha = { 0.01, 0.02, 0.03 };       
    beta  = { 0.0, 0.006, 0.01, 0.02 };
    temp  = { 296.0 };                                          
    delta = 0.00255;                                           
    twt    = 0.055556;     c      = 0.0;    tbeta = 0.444444;   
    oscE   = { 35.8,    0.48};                                 
    oscW   = { 0.166667, 0.333333 };                            
    dka    = 0.0;                            
    kappa  = { };                            

//4.5705275051479635E-009   4.5970339986086680E-009   4.6147049942491374E-009   4.6588824833503055E-009   9.1372603005644147E-009   9.1902512709537933E-009   9.2255785865584117E-009   9.3138968826866152E-009   1.3701367965741643E-008   1.3780825055741895E-008   1.3833796647506908E-008   1.3966226328467962E-008

    WHEN( "Energy grid is not specified (uniform grid assumed)" ){
      rho   = { 0.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 0.0 };
      auto out = leapr( nphon, awr, 
          ncold, aws, lat, alpha, beta, 
          temp, delta, rho, twt, c, tbeta, oscE, 
          oscW, dka, kappa );
  
      auto ssm = std::get<2>(out);
      std::vector<double> ssmCorrect { };
      //checkSab( ssmCorrect, ssm );
    } // WHEN

    WHEN( "Energy grid is specified (uniform grid not assumed)" ){
      rho = { 0.0, 2.0, 2.0, 3.0, 3.0, 0.0 };
      std::vector<double> energyGrid { 0.0, 2.0, 5.0, 6.0, 8.0, 9.0 };
      auto out = leapr( nphon, awr, 
          ncold, aws, lat, alpha, beta, 
          temp, delta, rho, twt, c, tbeta, oscE, 
          oscW, dka, kappa, energyGrid );
  
      auto ssm = std::get<2>(out);
      std::vector<double> ssmCorrect { };
      //checkSab( ssmCorrect, ssm );

    } // WHEN

  } // GIVEN 

} // TEST CASE
