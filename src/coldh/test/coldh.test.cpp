#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "discre/discre_util/bfill.h"
#include "discre/discre_util/exts.h"
#include "coldh/coldh.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

//auto equal = [](auto x, auto y, double tol = 5e-5){return x == Approx(y).epsilon(tol);};

TEST_CASE( "coldh" ){
  std::vector<double> tempf, alpha, beta, ska, goodSymSab1(25), goodSymSab2(25);
  double temp, tev, tbeta, trans_weight, scaling, dka;
  int ncold, lat;
  bool free;

  tempf = {193093.99765}; 
  alpha = {0.1, 0.2, 0.4, 0.8, 1.6},
  beta  = {0.10, 0.15, 0.30, 0.60, 1.2}; 
  ska   = { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13 };

  temp = 200.0, tev = 1.72346606E-2, tbeta = 0.7, trans_weight = 0.3, scaling = 1.4679720;

  double tbart = 2.1908278997651038; // tempf[t]/temp
  ncold = 1, lat = 1;

  std::vector<double> sym_sab_1 = ranges::view::iota(1,26);
  std::vector<double> sym_sab_2 (25,0.0);


  GIVEN( "molecular translations are assumed to not be free" ){
    free = false;

    WHEN( "step size is moderate" ){
      dka = 0.2;

      goodSymSab1 = { 1.31298597, 2.61610678, 3.89448375, 5.18220084, 6.37976960, 7.39472709, 8.70742257, 9.93419554, 11.2012531, 12.2744188, 12.3366109, 13.7116189, 14.8544522, 16.1144774, 16.9761810, 14.9544688, 16.5017425, 17.5200167, 18.8153504, 19.3560515, 14.0181895, 15.8631250, 16.7233471, 18.1109384, 18.2909144 };

      goodSymSab2 = { 1.31298597, 2.07970231, 2.48840556, 2.12146489, 1.11262046, 7.39472709, 6.93072507, 6.35851602, 4.62138876, 2.15202264, 12.3366109, 10.8627057, 9.49252611, 6.69514723, 3.01827095, 14.9544688, 12.9372515, 11.1527435, 7.92220342, 3.53796519, 14.0181895, 12.1754199, 10.5638393, 7.82654458, 3.50708747 };


      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, trans_weight, tbeta, scaling, 
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2, tbart );
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN

    /*

    WHEN( "step size is small" ){
      dka = 0.0001;

      goodSymSab1 = { 1.3287711, 2.6267533, 3.9235595, 5.2043802, 6.5088037, 
        7.6149802, 8.8664244, 10.115396, 11.306691, 12.613772, 13.088162, 
        14.251586, 15.410194, 16.408195, 17.759561, 16.788298, 17.794618, 
        18.792099, 19.406237, 20.903111, 17.308492, 18.066429, 18.809744, 
        18.777265, 20.590641 };

      goodSymSab2 = { 1.3287711, 2.2621297, 2.9081236, 2.8633763, 1.9675888,
        7.6149802, 7.6447762, 7.5170871, 6.2697456, 3.8402507, 13.088162, 
        12.312841, 11.499779, 9.2237236, 5.4789194, 16.788298, 15.438241, 
        14.144248, 11.221266, 6.6252107, 17.308492, 15.813446, 14.417152, 
        11.526322, 6.8941688 };

      THEN( "output scattering law vectors are correct" ){
        coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
//        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
 //       REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
 
      AND_WHEN( "kappa values are smaller" ){
        ska = { 0.1, 0.2, 0.3, 0.5, 0.8, 1.13 };

        THEN( "output scattering law vectors are correct" ){
          coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
            alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
  //          REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
   //         REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
        } // THEN
      } // AND WHEN
    } // WHEN


    WHEN( "step size is large" ){
      dka = 10.0;
      ska = { 100.001, 5.002, 0.003 };

      goodSymSab1 = { 4.712662, 9.394536, 14.07523, 18.73994, 23.42826, 
        26.36307, 30.73920, 35.11286, 39.42884, 43.86060, 43.00731, 46.89065, 
        50.76919, 54.48711, 58.55840, 50.71420, 53.84089, 56.95874, 59.69325, 
        63.31049, 45.48792, 47.58773, 49.67292, 50.98232, 54.137579 };

      goodSymSab2 = { 4.712662, 8.087215, 10.42866, 10.29185, 7.063631, 
        26.36307, 26.47085, 26.03566, 21.70350, 13.25161, 43.00731, 40.40555, 
        37.69436, 30.12187, 17.76729, 50.71420, 46.46355, 42.41879, 33.33124, 
        19.39806, 45.48792, 41.22266, 37.28116, 29.20083, 16.998312 };

      THEN( "output scattering law vectors are correct" ){
        coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
//        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
 //       REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "molecular translations are assumed to be free" ){
    free = true; ska = { 1.1, 2.2, 3.3, 4.4 }; dka = 2.3; tempf[0] = 2e5;

    WHEN( "ortho hydrogen is chosen" ){
      ncold = 1; 

      goodSymSab1 = { 0.6906337, 0.7707476, 0.7716237, 0.6687087, 0.2834944, 
        0.4527882, 0.5086638, 0.5282524, 0.5302170, 0.4020220, 0.2736184, 
        0.3085437, 0.3264013, 0.3525404, 0.3571662, 0.1412989, 0.1596664, 
        0.1704781, 0.1909689, 0.2234552, 5.380649E-2, 6.086823E-2, 6.529124E-2, 
        7.447846E-2, 9.362654E-2 };
    
      goodSymSab2 = { 0.6906337, 0.6648041, 0.5741337, 0.3704882, 
        8.711346E-2,  0.4527882, 0.4388417, 0.3931755, 0.2936692, 0.1230125, 
        0.2736184, 0.2660678, 0.2427135, 0.1949027, 0.1090072, 0.1412989, 
        0.1376181, 0.1266453, 0.1053845, 6.801453E-2, 5.380649E-2, 5.244522E-2, 
        4.847128E-2, 4.104705E-2, 2.843402E-2 };

      THEN( "output scattering law vectors are correct" ){
        coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
  //      REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
   //     REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN

    WHEN( "para hydrogen is chosen" ){
      ncold = 2;

      goodSymSab1 = { 3.761848E-2, 4.619964E-2, 5.011692E-2, 5.386082E-2, 
        4.389428E-2, 3.985561E-2, 4.787329E-2, 5.246895E-2, 6.015803E-2, 
        6.454883E-2, 4.580030E-2, 5.319663E-2, 5.767247E-2, 6.629526E-2, 
        7.942194E-2, 4.741764E-2, 5.420193E-2, 5.841111E-2, 6.698852E-2, 
        8.363218E-2, 3.549853E-2, 4.035129E-2, 4.343504E-2, 4.996043E-2, 
        6.417559E-2 };

      goodSymSab2 = { 3.761848E-2, 3.572692E-2, 2.999425E-2, 1.959557E-2, 
        8.295755E-3, 3.985561E-2, 3.826508E-2, 3.363161E-2, 2.537094E-2, 
        1.394832E-2, 4.580030E-2, 4.435546E-2, 4.013064E-2, 3.232050E-2, 
        1.984000E-2, 4.741764E-2, 4.610451E-2, 4.226779E-2, 3.511375E-2, 
        2.316111E-2, 3.549853E-2, 3.457250E-2, 3.188555E-2, 2.692709E-2, 
        1.866112E-2 };
    
      THEN( "output scattering law vectors are correct" ){
        coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
//        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
 //       REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
 
    } // WHEN

    WHEN( "ortho deuterium is chosen" ){
      ncold = 3;

      goodSymSab1 = { 1.380304, 1.543119, 1.546022, 1.340770, 0.5632300, 
        0.9730467, 1.095199, 1.138230, 1.142990, 0.8615740, 0.6178966, 
        0.6978264, 0.7386705, 0.7982229, 0.8065778, 0.3149662, 0.3563194, 
        0.3806761, 0.4268510, 0.4999492, 0.1121509, 0.1269906, 0.1363143, 
        0.1557608, 0.1967059 };

      goodSymSab2 = { 1.380304, 1.328098, 1.145199, 0.7357275, 0.1696729, 
        0.9730467, 0.9426315, 0.8431951, 0.6272528, 0.2594906, 0.6178966, 
        0.6006209, 0.5472136, 0.4380633, 0.2429255, 0.3149662, 0.3066842, 
        0.2820068, 0.2342521, 0.1505688, 0.1121509, 0.1092965, 0.1009743, 
        8.546640E-2, 5.922167E-2 };
      
    
      THEN( "output scattering law vectors are correct" ){
        coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
//        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
 //       REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
 
    } // WHEN

    WHEN( "para deuterium is chosen" ){
      ncold = 4;

      goodSymSab1 = { 1.288188, 1.439972, 1.442868, 1.252397, 0.5294604, 
        0.9198252, 1.035278, 1.076075, 1.081155, 0.8172565, 0.5941633, 
        0.6710235, 0.7103489, 0.7678541, 0.7769734, 0.3101028, 0.3508164, 
        0.3748090, 0.4203312, 0.4926102, 0.1138499, 0.1289094, 0.1383714, 
        0.1581061, 0.1996595 };

      goodSymSab2 = { 1.288188, 1.239549, 1.069146, 0.6875359, 0.1594064, 
        0.9198252, 0.8911030, 0.7972289, 0.5934155, 0.2461694, 0.5941633, 
        0.5775611, 0.5262499, 0.4214230, 0.2340327, 0.3101028, 0.3019515, 
        0.2776671, 0.2306855, 0.1483746, 0.1138499, 0.1109532, 0.1025076, 
        8.676988E-2, 6.013533E-2 };
      
    
      coldh( temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2 );
//        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
 //       REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
 
    } // WHEN
*/
  } // GIVEN
} // TEST CASE
