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
    if (aInit == aFinal){
      REQUIRE( correct[i] == Approx(ssm(aInit,bInit+i,0)).epsilon(tol) ); 
    }
    else if (bInit == bFinal){
      REQUIRE( correct[i] == Approx(ssm(aInit+i,bInit,0)).epsilon(tol) ); 
    }
    else {
      REQUIRE( correct[i] == Approx(ssm(aInit+i,bInit+i,0)).epsilon(tol) ); 
    }
  } 

}




TEST_CASE( "leapr" ){
  int nout, ntempr, iprint, nphon, mat, npr, iel, ncold, nss, nalpha, nbeta, 
      lat, ni, nd, nka;
  //int mss;
  double sps, b7, za, awr, spr, aws, delta, twt, c, tbeta, dka;
  std::vector<double> alpha, beta, temp, rho, oscE, oscW, 
    kappa;
  std::string title;

  GIVEN( "simple H in H20 input" ) {
    nout   = 24;                                                 // Card 1
    title  = "h in h20, shortened endf model";                   // Card 2
    ntempr = 1;       iprint = 1;      nphon = 100;              // Card 3
    mat    = 101;     za     = 1001.0;                           // Card 4
    awr    = 0.99917; spr    = 20.449; npr   = 2;   iel = 0;   ncold = 0; // Card 5
    nss    = 1;       b7     = 1.;     aws   = 1.1; sps = 3.8883; //mss = 1; 
                                                                     // Card 6
    nalpha = 5;       nbeta  = 7;      lat   = 1;                // Card 7
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406 };        // Card 8
    beta  = { 0.0, 0.006375, 0.01275, 0.0255, 0.03825, 0.051, 0.06575 }; // Card 9
    temp  = { 296.0 };                                          // Card 10 
    delta = 0.00255;     ni = 67;                           // Card 11
    rho   = { 0.0000, 0.0005, 0.0010, 0.0020, 0.0035, 0.0050, 0.0075, 0.0100, 
      0.0130, 0.0165, 0.0200, 0.0245, 0.0290, 0.0340, 0.0395, 0.0450, 0.0506, 
      0.0562, 0.0622, 0.0686, 0.0750, 0.0830, 0.0910, 0.0990, 0.1070, 0.1150, 
      0.1197, 0.1214, 0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.0954, 0.0912, 
      0.0871, 0.0839, 0.0807, 0.0779, 0.0757, 0.0735, 0.0716, 0.0697, 0.0680, 
      0.0665, 0.0650, 0.0634, 0.0618, 0.0602, 0.0586, 0.0571, 0.0558, 0.0546, 
      0.0535, 0.0525, 0.0515, 0.0504, 0.0493, 0.0482, 0.0470, 0.0459, 0.0447, 
      0.0436, 0.0428, 0.0424, 0.0420, 0.0000 };
                                                                     // Card 12
    twt    = 0.055556;     c      = 0.0;    tbeta = 0.444444;    // Card 13
    nd     = 2;                                                  // Card 14
    oscE   = { 35.8,    0.48};                                  // Card 15
    oscW   = { 0.166667, 0.333333 };                             // Card 16
    nka    = 4;       dka    = 0.01;                             // Card 17
    kappa  = { 0.1, 0.2, 0.4, 0.7 };                             // Card 18

    auto out = leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, 
        spr, npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, 
        temp, delta, ni, rho, twt, c, tbeta, nd, oscE, 
        oscW, nka, dka, kappa );

    auto ssm = std::get<2>(out);
    std::vector<double> ssmCorrect { 11.9380466, 11.7358569, 11.1552444,  
      9.04435365, 6.36253466, 3.86288069, 1.81640727, 9.77386851, 9.67827226, 
      9.35649187, 8.13748192, 6.42985286, 4.61578595, 2.78373847, 7.52086114, 
      7.48096532, 7.33830900, 6.77191769, 5.91039840, 4.85809341, 3.60668537, 
      6.55904026, 6.53116851, 6.44104253, 6.08124304, 5.47109680, 4.72905752, 
      3.76787280, 5.28338485, 5.26835188, 5.23245837, 5.03910024, 4.71891225,  
      4.29716228, 3.7169103340 };

    checkSab( ssmCorrect, ssm );
  } // GIVEN 

  GIVEN( "H in H2O input (TEST 09)" ) {

  nout   = 20;                                                      // Card 1
  title  = "H IN H2O, ENDF MODEL";                                  // Card 2
  ntempr = 1;       iprint = 1;     nphon = 100;                    // Card 3
  mat    = 1;       za     = 1001.0;                                // Card 4
  awr    = 0.99917; spr    = 20.478;  npr   = 2;  iel = 0;   ncold = 0;
  nss    = 1;       b7     = 1.0;   aws = 15.85316;   sps = 3.8883; // Card 6
  nalpha = 65;      nbeta  = 75;    lat   = 1;                      // Card 7
  alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 0.100812, 0.151218,
    0.201624, 0.25203, 0.302436, 0.352842, 0.403248, 0.453654, 0.50406, 
    0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
    1.07213, 1.17708, 1.29211, 1.41822, 1.55633, 1.70775, 1.87379, 2.05566, 
    2.25506, 2.47352, 2.71295, 2.97546, 3.26308, 3.57832, 3.92390, 4.30266, 
    4.71770, 5.17256, 5.67118, 6.21758, 6.81650, 7.47289, 8.19228, 8.98073, 
    9.84489, 10.7919, 11.8303, 12.9674, 14.2145, 15.5815, 17.0796, 18.7208, 
    20.5203, 22.4922, 24.6526, 27.0216, 29.6175, 32.4625, 35.5816, 38.9991, 
    42.7453, 46.8503, 50.0 };  // Card 8
  beta = { 0.00000, 0.006375, 0.01275, 0.02550, 0.03825, 0.05100, 0.06575, 
    0.0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 
    0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.04844, 
    1.12909, 1.20974, 1.29039, 1.37104, 1.45169, 1.53234, 1.61299, 1.69364, 
    1.77429, 1.85494, 1.93559, 2.01624, 2.09689, 2.17754, 2.25819, 2.33884, 
    2.41949, 2.50014, 2.58079, 2.66950, 2.76709, 2.87445, 2.99250, 3.12235, 
    3.26530, 3.42247, 3.59536, 3.78549, 3.99467, 4.22473, 4.47787, 4.75631, 
    5.06258, 5.39939, 5.76997, 6.17766, 6.62607, 7.11924, 7.66181, 8.25862, 
    8.91511, 9.63722, 10.4320, 11.3051, 12.2668, 13.3243, 14.4867, 15.0 }; // Card 9
  temp   = { 296.0 };                                             // Card 10 
  delta  = 0.00255;    ni = 67;                                   // Card 11
  rho = { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 
    0.013, 0.0165, 0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 
    0.0622, 0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
    0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 0.0839, 
    0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 0.06652, 0.065, 
    0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 0.05462, 0.0535, 0.0525, 
    0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 0.0459, 0.04478, 0.04366, 
    0.04288, 0.04244, 0.042, 0.0 };           // Card 12
  twt = 0.055556; c = 0.0; tbeta = 0.444444;  // Card 13
  nd = 2;                                     // Card 14
  oscE = { 0.205, 0.48 };                     // Card 15
  oscW = { 0.166667, 0.333333 };              // Card 16
  nka = 0; dka = 0.0;                         // Card 17
  kappa    = { };                             // Card 18

  auto out = leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, 
      spr, npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, 
      temp, delta, ni, rho, twt, c, tbeta, nd, oscE, 
      oscW, nka, dka, kappa );

  double lambda_s = std::get<0>(out);
  double t_eff = std::get<1>(out);
  auto ssm = std::get<2>(out);

  std::vector<double> ssm_alpha_0_49_beta_0_49 { 11.935614, 9.6753951, 7.3348118, 6.0777375, 
    4.7152145, 3.7635737, 3.1352778, 2.5435978, 1.9451342, 1.5299730, 0.9617595, 
    0.5859237, 0.3524469, 0.2124005, 0.1300019, 8.2531345e-2, 5.6886704e-2, 
    4.3458881e-2, 3.6822232e-2, 3.4067983e-2, 3.3632585e-2, 3.4628352e-2, 
    3.6546693e-2, 3.9131910e-2, 4.2191691e-2, 4.5634099e-2, 4.9388018e-2, 
    5.3368441e-2, 5.7527799e-2, 6.1792375e-2, 6.6065929e-2, 7.0233475e-2, 
    7.4218351e-2, 7.7872169e-2, 8.1100605e-2, 8.3812876e-2, 8.5914606e-2, 
    8.7365248e-2, 8.8101263e-2, 8.8115627e-2, 8.7348737e-2, 8.5804720e-2, 
    8.3492917e-2, 8.0470802e-2, 7.6789248e-2, 7.2560911e-2, 6.7879552e-2, 
    6.2878926e-2, 5.7644932e-2, 5.2325611e-2, 4.6996838e-2, 4.1750932e-2, 
    3.6686082e-2, 3.1863917e-2, 2.7336267e-2, 2.3145924e-2, 1.9324611e-2, 
    1.5889523e-2, 1.2853166e-2, 1.0214210e-2, 7.9612973e-3, 6.0777445e-3, 
    4.5355538e-3, 3.3038095e-3, 2.7361841e-3 };
  std::vector<double> ssm_alpha_29_beta_0_39 { 0.5075435, 0.5087149, 0.5099015, 
    0.5120805, 0.5143626, 0.5168548, 0.5195800, 0.5205411, 0.5222588, 0.5201505, 
    0.5066769, 0.4814219, 0.4461412, 0.4027074, 0.3550043, 0.3064586, 0.2598735, 
    0.2170864, 0.1793159, 0.1473404, 0.1215918, 0.1017288, 8.700412e-2, 
    7.656394e-2, 6.950085e-2, 6.510340e-2, 6.265127e-2, 6.156786e-2, 6.138864e-2, 
    6.179237e-2, 6.254617e-2, 6.347645e-2, 6.448598e-2, 6.546971e-2, 6.638564e-2, 
    6.717356e-2, 6.778630e-2, 6.821271e-2, 6.839302e-2, 6.832765e-2 };
  std::vector<double> ssm_alpha_0_39_beta_40 { 5.57415136e-4, 8.28335298e-4, 
    1.38765192e-3, 1.81319262e-3, 2.75597276e-3, 4.10379550e-3, 5.43260358e-3, 
    8.03115066e-3, 1.05535283e-2, 1.30041386e-2, 1.53847306e-2, 1.77000686e-2, 
    1.99510266e-2, 2.21390179e-2, 2.42686023e-2, 2.63383940e-2, 2.85428116e-2, 
    3.08854175e-2, 3.33695419e-2, 3.59893717e-2, 3.87456221e-2, 4.16348357e-2, 
    4.46483735e-2, 4.77815183e-2, 5.10161200e-2, 5.43383922e-2, 5.77140508e-2, 
    6.11324854e-2, 6.45536433e-2, 6.79373750e-2, 7.12310538e-2, 7.43908691e-2, 
    7.73717267e-2, 8.01064016e-2, 8.25233157e-2, 8.45724972e-2, 8.61959064e-2, 
    8.73502727e-2, 8.79674594e-2, 8.79771268e-2 };
   
  checkPartofSab(  0, 49,  0, 49, ssm_alpha_0_49_beta_0_49, ssm, 5e-5 );
  checkPartofSab( 29, 29,  0, 39, ssm_alpha_29_beta_0_39,   ssm, 1e-6 );
  checkPartofSab(  0, 39, 40, 40, ssm_alpha_0_39_beta_40,   ssm, 1e-6 );
 
  REQUIRE( 0.23520650 == Approx(lambda_s).epsilon(1e-6) );        
  REQUIRE( 1.93448465 == Approx(t_eff).epsilon(1e-6) );        

  b7  += 1.0;
  sps += 1.0;

 
  } // GIVEN 


  GIVEN( "Para hydrogen at 20K (TEST 22)" ) {

    nout   = 20;                                                    // Card 1
    title  = "para hydrogen at 20k";                                // Card 2
    ntempr = 1;       iprint = 2;     nphon = 100;                  // Card 3
    mat    = 2;       za     = 1001.0;  
    // int isabt = 0; int ilog = 0; 
    double smin = 2.0e-38;                                           // Card 4
    awr    = 0.99917; spr   = 20.478; npr = 2;   iel = 0; ncold = 2; // Card 5
    nss    = 0;       b7    = 0.0;    aws = 0.0; sps = 0.0;          // Card 6
    nalpha = 59;      nbeta = 105;    lat = 0;                       // Card 7
    std::vector<double> alpha = { 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 
      0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 
      14, 16, 18, 20, 25, 30, 40, 50, 60, 80, 100, 120, 140, 160, 180, 200,
      220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 
      500, 520, 540, 560, 580, 600 };                                // Card 8
    std::vector<double> beta = { 0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 
      0.4, 0.5, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 5.7, 
      5.8, 5.9, 6, 6.25, 6.5, 7, 7.5, 8, 8.4, 8.5, 8.6, 8.7027, 8.8, 
      8.9, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 
      15, 15.5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 25.7, 26.1, 
      26.5, 27, 28, 29, 30, 35, 37.5, 40, 45.5, 50, 52.2, 55, 60, 65, 
      70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 
      140, 145, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 
      260, 280, 300 }; // Card 9
    temp   = { 20.0 };                                             // Card 10 
    delta  = 0.00025;    ni = 48;                                  // Card 11
    std::vector<double> rho = { 0.0, 0.01563, 0.0625, 0.141, 0.25, 
      0.391, 0.5625, 0.766, 1., 1.266, 1.5625, 1.89, 2.25, 2.64, 3.0625, 
      3.52, 4., 4.6, 5.5, 7.0, 8.5, 9.2, 9.5, 9.4, 9.2, 8.9, 8.5, 8.0, 
      7.5, 7.05, 6.7, 6.4, 6.2, 6.1, 6.2, 6.45, 6.7, 6.95, 7.1, 6.55, 5.5, 
      3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };      // Card 12
    twt = 0.025; c = 40.0; tbeta = 0.475;       // Card 13
    nd = 0;                                     // Card 14
    oscE = { };                                 // Card 15
    oscW = { };                                 // Card 16
    nka = 301; dka = 0.05;                      // Card 17
    kappa = { 9.29815e-2, 9.29815e-2, 9.31478e-2, 9.43837e-2, 9.5458e-2, 
      9.68871e-2, 9.86644e-2, 1.00813e-1, 1.03358e-1, 1.06334e-1, 1.09777e-1, 
      1.1373e-1, 1.18249e-1, 1.23395e-1, 1.29241e-1, 1.35875e-1, 1.43399e-1, 
      1.51932e-1, 1.61618e-1, 1.72622e-1, 1.85143e-1, 1.99416e-1, 2.15717e-1, 
      2.34376e-1, 2.55779e-1, 2.80385e-1, 3.08731e-1, 3.41443e-1, 3.7924e-1, 
      4.22934e-1, 4.73417e-1, 5.3162e-1, 5.98441e-1, 6.74614e-1, 7.9867e-1, 
      1.118, 1.4087, 1.6633, 1.8799, 2.0554, 2.1889, 2.2803, 2.3306, 2.3418, 
      2.3164, 2.2581, 2.1706, 2.0585, 1.9264, 1.779, 1.6213, 1.4581, 1.2938, 
      1.1327, 9.7879e-1, 8.3545e-1, 7.0565e-1, 5.918e-1, 4.9576e-1, 4.1882e-1, 
      3.6172e-1, 3.2464e-1, 3.0727e-1, 3.0884e-1, 3.2812e-1, 3.6355e-1, 
      4.1325e-1, 4.751e-1, 5.4681e-1, 6.2601e-1, 7.1023e-1, 7.9708e-1, 8.8421e-1, 
      9.6942e-1, 1.0507, 1.1261, 1.1943, 1.2538, 1.3037, 1.3432, 1.3719, 1.3897, 
      1.3968, 1.3935, 1.3805, 1.3586, 1.3289, 1.2924, 1.2505, 1.2045, 1.1557, 
      1.1055, 1.0551, 1.0059, 9.5899e-1, 9.1541e-1, 8.7608e-1, 8.4174e-1, 8.1298e-1, 
      7.9025e-1, 7.738e-1, 7.6372e-1, 7.5995e-1, 7.6226e-1, 7.7032e-1, 7.8364e-1, 
      8.0164e-1, 8.2367e-1, 8.4898e-1, 8.7681e-1, 9.0637e-1, 9.3684e-1, 9.6745e-1, 
      9.9744e-1, 1.0261, 1.0529, 1.0771, 1.0984, 1.1163, 1.1306, 1.141, 1.1476, 
      1.1504, 1.1493, 1.1447, 1.1369, 1.126, 1.1126, 1.0971, 1.0799, 1.0616, 
      1.0426, 1.0233, 1.0044, 9.8614e-1, 9.6901e-1, 9.5336e-1, 9.3949e-1, 9.2766e-1, 
      9.1806e-1, 9.1082e-1, 9.0602e-1, 9.0365e-1, 9.0367e-1, 9.0598e-1, 9.1043e-1, 
      9.1682e-1, 9.2491e-1, 9.3444e-1, 9.4513e-1, 9.5667e-1, 9.6876e-1, 9.8109e-1, 
      9.9335e-1, 1.0053, 1.0166, 1.027, 1.0364, 1.0445, 1.0513, 1.0566, 1.0603, 
      1.0625, 1.0631, 1.0622, 1.0599, 1.0563, 1.0515, 1.0457, 1.0391, 1.0318, 1.0241, 
      1.0161, 1.008, 1.0001, 9.925e-1, 9.8538e-1, 9.7887e-1, 9.7311e-1, 9.682e-1, 
      9.6422e-1, 9.6122e-1, 9.5924e-1, 9.5828e-1, 9.5832e-1, 9.5932e-1, 9.6122e-1, 
      9.6394e-1, 9.6739e-1, 9.7145e-1, 9.7602e-1, 9.8096e-1, 9.8614e-1, 9.9145e-1, 
      9.9674e-1, 1.0019, 1.0068, 1.0114, 1.0155, 1.0192, 1.0222, 1.0246, 1.0264, 
      1.0275, 1.0279, 1.0277, 1.0268, 1.0254, 1.0234, 1.021, 1.0182, 1.0151, 1.0117, 
      1.0082, 1.0047, 1.0012, 9.9774e-1, 9.945e-1, 9.9151e-1, 9.8882e-1, 9.8648e-1, 
      9.8454e-1, 9.8302e-1, 9.8193e-1, 9.813e-1, 9.8111e-1, 9.8135e-1, 9.820e-1, 
      9.8303e-1, 9.8441e-1, 9.8608e-1, 9.8801e-1, 9.9013e-1, 9.924e-1, 9.9475e-1, 
      9.9713e-1, 9.9949e-1, 1.0018, 1.0039, 1.0059, 1.0077, 1.0092, 1.0105, 1.0115, 
      1.0122, 1.0126, 1.0127, 1.0125, 1.012, 1.0113, 1.0104, 1.0092, 1.0079, 1.0065, 
      1.0049, 1.0034, 1.0017, 1.0001, 9.986e-1, 9.9715e-1, 9.958e-1, 9.946e-1, 
      9.9355e-1, 9.9269e-1, 9.9201e-1, 9.9154e-1, 9.9126e-1, 9.9119e-1, 9.9132e-1, 
      9.9163e-1, 9.9211e-1, 9.9274e-1, 9.9351e-1, 9.944e-1, 9.9537e-1, 9.9641e-1, 
      9.9749e-1, 9.9859e-1, 9.9967e-1, 1.0007, 1.0017, 1.0026, 1.0035, 1.0042, 
      1.0048, 1.0053, 1.0056, 1.0058, 1.0058, 1.0058, 1.0056, 1.0053, 1.0048, 1.0043, 
      1.0037, 1.0031, 1.0024, 1.0016, 1.0009, 1.0001, 9.994e-1 };      // Card 18

    auto out = leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, spr, 
      npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, temp, delta, 
      ni, rho, twt, c, tbeta, nd, oscE, oscW, nka, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);
    REQUIRE( 0.18256462793803138 == Approx(lambda_s).epsilon(1e-6) );        

    std::vector<double> ssm_alpha_0_30_beta_0_30 { 0.61800579, 4.0317687e-2, 
      1.8361451e-3, 2.4153060e-4, 6.9097368e-5, 7.4382647e-5, 2.2143760e-5, 
      1.1852629e-5, 1.2697481e-5, 2.0489166e-5, 4.1755430e-5, 7.6980980e-5, 
      1.4935133e-4, 3.5059352e-4, 7.3923986e-4, 1.7171400e-3, 4.4357765e-3, 
      6.2911089e-3, 6.6963705e-3, 6.4565419e-3, 5.7453644e-3, 4.3886776e-3, 
      4.2063259e-3, 4.7262951e-3, 5.6954504e-3, 6.3167771e-3, 7.8156250e-3, 
      6.6440418e-3, 5.0127895e-3, 1.9578449e-3, 6.2127409e-4 };
    std::vector<double> ssm_alpha_30_58_beta_58 { 3.3474674e-2, 1.7936553e-2, 
      3.3308960e-3, 4.1497349e-4, 3.7854146e-5, 2.5638317e-6, 1.5570105e-7, 
      2.4281266e-8, 6.9117027e-9, 1.6950933e-9, 3.4247961e-10, 5.9755048e-11, 
      9.3828990e-12, 1.3583384e-12, 1.8392439e-13, 2.3513069e-14, 2.8553552e-15, 
      3.3078253e-16, 3.6727172e-17, 3.9361536e-18, 4.1320631e-19, 4.3767425e-20, 
      4.9157113e-21, 6.1841146e-22, 8.8819050e-23, 1.4063636e-23, 2.3241601e-24, 
      3.8497620e-25, 6.2600002e-26 };
    std::vector<double> ssm_alpha_30_beta_70_104 { 1.0700481e-2, 7.9260893e-3, 
      5.5082711e-3, 3.6493970e-3, 4.0224294e-3, 5.0148042e-3, 5.2654162e-3, 
      4.4472258e-3, 3.0007869e-3, 1.7773242e-3, 8.7339023e-4, 3.9012523e-4, 
      1.6011654e-4, 6.2585891e-5, 2.5439285e-5, 1.0846789e-5, 4.8935073e-6, 
      3.2775460e-6, 7.2170171e-6, 1.7372343e-5, 2.9119281e-5, 3.5241659e-5, 
      2.3687997e-5, 7.7156439e-6, 1.4480755e-6, 1.7845188e-7, 4.5460306e-8, 
      2.5256187e-8, 8.6495516e-9, 1.6957101e-9, 1.3974824e-9, 2.0532145e-8, 
      4.7822041e-8, 1.2107986e-8, 2.8016681e-10 };

    checkPartofSab(  0, 30,  0, 30,  ssm_alpha_0_30_beta_0_30, ssm, 1e-6 );
    checkPartofSab( 30, 58, 58, 58,  ssm_alpha_30_58_beta_58,  ssm, 1e-6 );
    checkPartofSab( 30, 30, 70, 104, ssm_alpha_30_beta_70_104, ssm, 1e-6 );
 

  } // GIVEN 



  GIVEN( "Beryllium Metal (pg 25 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
         "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
         "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {

    nout   = 20;                                                    // Card 1
    title  = "BE METAL, ENDF MODEL (EXTENDED)";                                // Card 2
    ntempr = 8;       iprint = 1;     nphon = 150;                  // Card 3
    mat    = 26;       za     = 126.0;  
    // int isabt = 0; int ilog = 0; 
    // double smin = 2.0e-38;                                       // Card 4
    awr    = 8.93478; spr   = 6.15; npr = 1;   iel = 2; ncold = 0; // Card 5
    nss    = 0;       b7    = 0.0;    aws = 0.0; sps = 0.0;          // Card 6
    nalpha = 90;      nbeta = 103;    lat = 1;                       // Card 7
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
              70.0, 75.0, 80.0 };                                // Card 8
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
            46.0, 48.0, 50.0, 52.5, 55.0, 57.5, 60.0, 65.0, 70.0, 75.0, 80.0 }; // Card 9
    temp   = { 296.0 };                                             // Card 10 
    delta  = 0.0015361;    ni = 55;                                  // Card 11
    rho = { 0.0, 0.0062, 0.025, 0.06, 0.115, 0.17, 0.235, 0.308, 0.39, 0.475, 
            0.595, 0.775, 1.052, 1.1482, 1.6233, 1.5653, 1.7423, 2.4289, 2.8072, 
            3.2863, 3.7577, 4.4397, 5.4924, 6.3315, 7.6421, 9.5339, 12.1016, 
            15.0553, 22.0154, 26.6382, 29.3387, 32.6036, 35.777, 37.6536, 
            40.3845, 37.1196, 32.39, 29.8269, 21.2831, 18.9488, 13.0186, 
            7.33848, 12.5075, 18.9031, 23.5564, 25.6923, 25.5855, 43.0544, 
            17.8656, 6.81975, 3.98963, 2.99184, 1.6843, 0.65772, 0.0 };  // Card 12
    twt = 0.0; c = 0.0; tbeta = 1.0;       // Card 13
    nd = 0;                                     // Card 14
    oscE = { };                                 // Card 15
    oscW = { };                                 // Card 16
    nka = 0; dka = 0.0;                      // Card 17
    kappa = { };      // Card 18

    auto out = leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, spr, 
      npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, temp, delta, 
      ni, rho, twt, c, tbeta, nd, oscE, oscW, nka, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    std::vector<double> ssm_alpha_0_20_beta_0_20 { 4.3597304e-4, 7.2685738e-4, 
      1.3786344e-3, 1.7585901e-3, 2.7464605e-3, 5.1608578e-3, 7.2133356e-3, 
      1.3417112e-2, 2.4827829e-2, 4.0356860e-2, 7.5381076e-2, 0.15624046, 
      0.25237716, 0.28748541, 0.30712752, 0.31311066, 0.28262440, 0.23217740, 
      0.18987750, 0.14735439, 0.11723088 };
    std::vector<double> ssm_alpha_89_beta_40_79 { 5.2517889e-9, 5.7583813e-9, 
      6.3410553e-9, 6.6687455e-9, 7.0125591e-9, 7.3915164e-9, 7.7904716e-9, 
      8.6956484e-9, 9.7542324e-9, 1.0997335e-8, 1.2464195e-8, 1.4203760e-8, 
      1.6279761e-8, 1.8769975e-8, 2.1772302e-8, 2.5414555e-8, 2.9864445e-8, 
      3.5327928e-8, 4.2085142e-8, 5.0506875e-8, 6.1032434e-8, 7.4341341e-8, 
      9.1266528e-8, 1.1292724e-7, 1.4091817e-7, 1.7735439e-7, 2.2503799e-7, 
      2.8810279e-7, 3.7217548e-7, 4.8506308e-7 };
    std::vector<double> ssm_alpha_30_70_beta_100 { 5.7623971e-22, 2.4059638e-21, 
      9.5534907e-21, 3.6214510e-20, 1.3139209e-19, 4.5780655e-19, 1.5349468e-18, 
      4.9643815e-18, 1.5521833e-17, 4.7009318e-17, 1.3806356e-16, 3.9398451e-16, 
      1.0946042e-15, 2.9577671e-15, 7.8061625e-15, 2.0078724e-14, 5.0489319e-14, 
      1.2414697e-13, 5.1105965e-13, 7.0492737e-13, 1.6315236e-12, 3.6978684e-12, 
      8.2270704e-12, 1.7964868e-11, 3.8573604e-11, 8.1330158e-11, 1.6887905e-10, 
      3.4443213e-10, 6.9185387e-10, 1.3200414e-9, 2.9481028e-9, 8.6594211e-9, 
      2.3585610e-8, 5.9991042e-8, 1.4337411e-7, 3.2368869e-7, 6.9358890e-7, 
      1.4164513e-6, 2.7671365e-6, 5.1882404e-6, 1.6310944e-5 };

    checkPartofSab(  0, 20,   0,  20, ssm_alpha_0_20_beta_0_20, ssm, 1e-6 );
    checkPartofSab( 30, 70, 100, 100, ssm_alpha_30_70_beta_100, ssm, 1e-6 );
    checkPartofSab( 89, 89,  40,  79, ssm_alpha_89_beta_40_79,  ssm, 1e-6 );
 
    REQUIRE( 0.73178450430084030 == Approx(lambda_s).epsilon(1e-6) );        

  } // GIVEN 

  GIVEN( "Beryllium Oxide (pg 39 of NEW THERMAL NEUTRON SCATTERING FILES FOR "
         "ENDF/B-VI RELEASE 2 by R.E. MacFarlane, 1994, available at "
         "https://t2.lanl.gov/nis/publications/thermal.pdf" ) {

    nout   = 20;                                                    // Card 1
    title  = "BEO, ENDF MODEL (EXTENDED)";                                // Card 2
    ntempr = 8;       iprint = 1;     nphon = 100;                  // Card 3
    mat    = 27;       za     = 127.0;  
    // int isabt = 0; int ilog = 0; 
    // double smin = 2.0e-38;                                       // Card 4
    awr    = 8.93478; spr   = 6.15; npr = 1;   iel = 3; ncold = 0; // Card 5
    nss    = 1;       b7    = 0.0;    aws = 15.858; sps = 3.7481; int mss = 1;          // Card 6
    nalpha = 90;      nbeta = 117;    lat = 1;                       // Card 7
    alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.0504, 0.0756, 0.1008, 0.15, 
              0.252, 0.33, 0.504, 0.756, 1.008, 1.260, 1.512, 1.764, 2.016, 
              2.268, 2.520, 2.772, 3.024, 3.282, 3.544, 3.813, 4.087, 4.366, 
              4.652, 4.943, 5.241, 5.545, 5.855, 6.172, 6.495, 6.825, 7.162, 
              7.507, 7.858, 8.217, 8.583, 8.957, 9.339, 9.729, 10.13, 10.53, 
              10.95, 11.37, 11.81, 12.25, 12.69, 13.16, 13.63, 14.11, 14.60, 
              15.10, 15.61, 16.13, 16.66, 17.21, 17.76, 18.3, 19.0, 20.0, 21.0, 
              22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 30.0, 32.0, 34.0, 36.0, 
              38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 52.5, 55.0, 57.5, 60.0, 
              62.5, 65.0, 70.0, 75.0, 80.0 };                                // Card 8
    beta = { 0.0, 0.1513, 0.3025, 0.4537, 0.6049, 0.7561, 0.9073, 1.059, 1.210, 
             1.361, 1.512, 1.5875, 1.663, 1.7390, 1.815, 1.8905, 1.966, 2.0415, 
             2.117, 2.1925, 2.268, 2.3435, 2.419, 2.4950, 2.571, 2.6465, 2.722, 
             2.7975, 2.873, 2.9485, 3.024, 3.1000, 3.176, 3.2515, 3.327, 3.4025, 
             3.478, 3.5850, 3.629, 3.7045, 3.780, 3.8560, 3.932, 4.0075, 4.083, 
             4.1620, 4.241, 4.3245, 4.408, 4.4955, 4.583, 4.6745, 4.766, 4.8620, 
             4.958, 5.0585, 5.159, 5.371, 5.592, 5.825, 6.069, 6.325, 6.593, 
             6.875, 7.170, 7.480, 7.805, 8.146, 8.504, 8.879, 9.273, 9.686, 
             10.12, 10.57, 11.05, 11.55, 12.07, 12.62, 13.20, 13.81, 14.44, 
             15.11, 15.81, 16.54, 17.31, 18.12, 18.96, 19.85, 20.78, 21.76, 
             22.78, 23.86, 24.99, 26.17, 27.41, 28.71, 30.08, 31.51, 33.01, 
             34.59, 36.24, 37.98, 39.80, 42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 
             54.0, 56.0, 58.0, 60.0, 65.0, 70.0, 75.0, 80.0 }; // Card 9
    temp   = { 296.0 };                                             // Card 10 
    delta  = 0.0016518;    ni = 84;                                  // Card 11
    rho = { 0.0, 0.3, 0.7, 0.9, 1.0, 1.2, 1.6, 2.0, 2.2, 3.0, 3.5, 4.5, 5.5, 
            6.8, 8.0, 9.2, 10.9, 12.9, 15.5, 18.6, 22.0, 26.0, 30.5, 35.0, 39.0, 
            40.0, 34.0, 28.0, 26.0, 24.4, 23.0, 21.3, 19.8, 17.0, 14.1, 12.0, 
            10.0, 9.0, 9.0, 8.5, 7.5, 6.0, 4.6, 3.1, 1.6, 0.5, 0.0, 0.0, 4.0, 
            15.0, 38.0, 52.0, 70.0, 105.0, 165.0, 230.0, 200.0, 170.0, 145.0, 
            136.0, 134.0, 112.0, 96.0, 89.0, 84.0, 75.0, 87.0, 81.0, 66.0, 59.0, 
            68.0, 105.0, 95.0, 97.0, 135.0, 163.0, 130.0, 111.0, 92.0, 67.0, 45.0, 
            19.0, 7.0, 0.0 };  // Card 12
    twt = 0.0; c = 0.0; tbeta = 1.0;       // Card 13
    nd = 0;                                     // Card 14
    oscE = { };                                 // Card 15
    oscW = { };                                 // Card 16
    nka = 0; dka = 0.0;                      // Card 17
    kappa = { };      // Card 18

    auto out = leapr( nout, title, ntempr, iprint, nphon, mat, za, awr, spr, 
      npr, iel, ncold, nss, aws, nalpha, nbeta, lat, alpha, beta, temp, delta, 
      ni, rho, twt, c, tbeta, nd, oscE, oscW, nka, dka, kappa );
  
    double lambda_s = std::get<0>(out), 
           t_eff    = std::get<1>(out);
    auto ssm = std::get<2>(out);

    std::vector<double> ssm_alpha_0_20_beta_0_20 { 2.7999522e-3, 2.2499360e-3, 
      1.4314522e-3, 1.5455335e-3, 2.2373451e-3, 3.6819532e-3, 5.5977873e-3, 
      9.2097916e-3, 1.8331881e-2, 2.8924029e-2, 4.9722500e-2, 6.9801132e-2, 
      7.5965863e-2, 6.8475570e-2, 6.4090814e-2, 6.2120218e-2, 6.0175214e-2, 
      5.7504955e-2, 5.4121428e-2, 4.9372518e-2, 4.5836374e-2 };
    std::vector<double> ssm_alpha_89_beta_90_116 { 4.2869564e-5, 5.7014993e-5, 
      7.6100497e-5, 1.0183761e-4, 1.3680959e-4, 1.8426649e-4, 2.4903698e-4, 
      3.3650597e-4, 4.5490864e-4, 6.1544286e-4, 8.3026609e-4, 1.1189152e-3, 
      1.5007918e-3, 2.0889495e-3, 2.7595450e-3, 3.5719780e-3, 4.5336106e-3, 
      5.6459639e-3, 6.9034170e-3, 8.2923692e-3, 9.7910059e-3, 1.1369440e-2, 
      1.2990680e-2, 1.6942691e-2, 2.0168093e-2, 2.2042809e-2, 2.2237244e-2 };
    std::vector<double> ssm_alpha_20_40_beta_116 { 7.8327127e-22, 3.4283860e-21, 
      1.3670939e-20, 5.0953897e-20, 1.7711576e-19, 5.7838000e-19, 1.7991308e-18, 
      5.3130914e-18, 1.5067278e-17, 4.1018293e-17, 1.0752489e-16, 2.7292713e-16, 
      6.7028902e-16, 1.6003523e-15, 3.7203584e-15, 8.4526482e-15, 1.8704134e-14, 
      4.0542011e-14, 8.5974181e-14, 1.7889490e-13, 3.655014e-13 };

    checkPartofSab(  0, 20,   0,  20, ssm_alpha_0_20_beta_0_20, ssm, 1e-6 );
    checkPartofSab( 20, 40, 116, 116, ssm_alpha_20_40_beta_116, ssm, 1e-6 );
    checkPartofSab( 89, 89,  90, 116, ssm_alpha_89_beta_90_116, ssm, 1e-6 );
 
    REQUIRE( 0.49335785747954103 == Approx(lambda_s).epsilon(1e-6) );        

  } // GIVEN 

} // TEST CASE
