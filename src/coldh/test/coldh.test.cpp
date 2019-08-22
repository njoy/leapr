#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "discre/discre_util/bfill.h"
#include "discre/discre_util/exts.h"
#include "coldh/coldh.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

//  std::cout << std::setprecision(15) << std::endl;
TEST_CASE( "coldh" ){

  std::vector<double> tempf, alpha, beta, ska, goodSymSab1(25), goodSymSab2(25);
  double temp, tev, continWeight, transWeight, scaling, dka;
  int ncold;
  bool free;

  tempf = {193093.99765}; 
  alpha = {0.1, 0.2, 0.4, 0.8, 1.6},
  beta  = {0.10, 0.15, 0.30, 0.60, 1.2}; 
  ska   = { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13 };

  temp = 200.0, tev = 1.72346606E-2, continWeight = 0.7, transWeight = 0.3, scaling = 1.4679720;
  for ( auto& a : alpha ){ a *= scaling; }
  for ( auto& b : beta  ){ b *= scaling; }

  double tbart = 2.1908278997651038; 
  

  dka = 0.2;
  std::vector<double> sym_sab_1 = ranges::view::iota(1,26);
  std::vector<double> sym_sab_2 (25,0.0);


  GIVEN( "Ortho hydrogen"){
    ncold = 1;

    WHEN( "molecular translations are assumed to not be free" ){
      free = false;
      goodSymSab1 = { 1.31298597, 2.61610678, 3.89448375, 5.18220084, 
      6.37976960, 7.39472709, 8.70742257, 9.93419554, 11.2012531, 12.2744188, 
      12.3366109, 13.7116189, 14.8544522, 16.1144774, 16.9761810, 14.9544688, 
      16.5017425, 17.5200167, 18.8153504, 19.3560515, 14.0181895, 15.8631250, 
      16.7233471, 18.1109384, 18.2909144 };

      goodSymSab2 = { 1.31298597, 2.07970231, 2.48840556, 2.12146489, 
      1.11262046, 7.39472709, 6.93072507, 6.35851602, 4.62138876, 2.15202264, 
      12.3366109, 10.8627057, 9.49252611, 6.69514723, 3.01827095, 14.9544688, 
      12.9372515, 11.1527435, 7.92220342, 3.53796519, 14.0181895, 12.1754199, 
      10.5638393, 7.82654458, 3.50708747 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN
    WHEN( "molecular translations are assumed to be free" ){
      free = true;
      goodSymSab1 = { 0.81319934, 0.92953946, 0.80962024, 0.37748227, 
      2.80179E-2, 0.54031344, 0.63051223, 0.62140318, 0.47593415, 0.12229667, 
      0.33102392, 0.39083244, 0.40988237, 0.40166856, 0.24921357, 0.17673224, 
      0.21016690, 0.22739072, 0.25166266, 0.24768537, 7.41382E-2, 8.85399E-2, 
      9.73032E-2, 0.11437214, 0.1419897 };

      goodSymSab2 = { 0.81319934, 0.74984840, 0.52789716, 0.16315535, 
      3.36375E-3, 0.54031344, 0.51003137, 0.40670334, 0.20395079, 2.09663E-2, 
      0.33102392, 0.31605061, 0.26796962, 0.17130854, 4.42535E-2, 0.17673224, 
      0.16963653, 0.14812323, 0.10666918, 4.40978E-2, 7.41382E-2, 7.13347E-2, 
      6.31585E-2, 4.81693E-2, 2.51138E-2 };
      
      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN
  } // GIVEN



  ska = { 0.1, 0.2, 0.3, 0.5, 0.8, 1.3, 2.1, 3.4, 5.5, 8.9, 14.4, 23.3, 37.7};
  dka = 0.5;


  GIVEN( "Para hydrogen"){
    ncold = 2;

    WHEN( "molecular translations are assumed to not be free" ){
      free = false;
      goodSymSab1 = { 0.1888934, 0.1668458, 0.1855598, 0.1452063, 0.4092368, 
      0.8484164, 0.8080851, 0.9374679, 1.1404978, 1.6914898, 2.4538030, 
      2.3868268, 2.8089883, 3.7222227, 4.8335209, 5.8832780, 5.7879925, 
      6.8452600, 9.3761663, 11.470245, 11.832878, 11.724655, 13.811441, 
      19.125832, 22.656718 };

      goodSymSab2 = { 0.1888934250, 0.1890980060, 0.1730912169, 0.1362329686, 
      2.3226989E-2, 0.8484164035, 0.8074767562, 0.7089605168, 0.5433885370, 
      0.1738668898, 2.4538030421, 2.3115056327, 2.0095825918, 1.5323268069, 
      0.5648625329, 5.8832780764, 5.5086753990, 4.7673249167, 3.6239109483, 
      1.4283935473, 11.832878126, 10.983505722, 9.4688868892, 7.1925924275, 
      2.9803987599 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN
    WHEN( "molecular translations are assumed to be free" ){
      free = true;
      goodSymSab1 = { 1.443219E-2, 2.968667E-2, 4.183377E-2, 6.001212E-2, 
      2.449593E-2, 3.029962E-2, 4.888518E-2, 6.180633E-2, 8.218185E-2, 
      6.525123E-2, 5.564104E-2, 7.404251E-2, 8.597772E-2, 0.107576001, 
      0.116909588, 7.754170E-2, 9.589133E-2, 0.107338492, 0.129571782, 
      0.159353849, 7.468211E-2, 9.041451E-2, 0.100414753, 0.121098379, 
      0.160501577 };

      goodSymSab2 = { 1.443219E-2, 1.234821E-2, 7.887972E-3, 5.747723E-3, 
      8.341121E-3, 3.029962E-2, 2.743585E-2, 2.081059E-2, 1.496107E-2, 
      1.132612E-2, 5.564104E-2, 5.238720E-2, 4.367380E-2, 3.079788E-2, 
      1.596106E-2, 7.754170E-2, 7.406808E-2, 6.416051E-2, 4.695455E-2, 
      2.294751E-2, 7.468211E-2, 7.170726E-2, 6.317298E-2, 4.800266E-2, 
      2.550662E-2 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "Ortho deuterium"){
    ncold = 3;
    WHEN( "molecular translations are assumed to not be free" ){
      free = false;
      AND_WHEN( "spacing (dka) is large" ){
        dka = 0.5;
        goodSymSab1 = { 1.04353549, 2.04283175, 3.04085742, 4.05944260, 
        5.05816298, 8.84191503, 10.3099118, 11.7696367, 13.2531960, 14.6990222, 
        27.5936205, 30.1836804, 32.7366103, 35.2069241, 37.9531925, 83.3216594, 
        89.4448115, 95.2689283, 99.3272795, 109.644553, 9.60975989, 10.6819098, 
        11.5852639, 11.7592124, 14.5146525 };

        goodSymSab2 = { 1.04353549, 1.64568291, 1.96673598, 1.67310377, 
        0.86654914, 8.84191503, 8.28251716, 7.58909537, 5.47516649, 2.51711392, 
        27.5936205, 24.2327449, 21.0875700, 14.5605371, 6.49983990, 83.3216594, 
        71.7845333, 61.3390811, 41.1194875, 18.7962242, 9.60975989, 8.58987526, 
        7.45990382, 4.82240089, 2.44514882 };

        THEN( "output scattering law vectors are correct" ){
          coldh( tev, ncold, transWeight + continWeight, 
            alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
          for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
            REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
            REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
          }
        } // THEN
      } // AND WHEN
      AND_WHEN( "spacing (dka) is small" ){
        dka = 0.05;
        goodSymSab1 ={ 1.6500137, 3.2545702, 4.8577548, 6.4805106, 8.0855785, 
        8.7834631, 10.241709, 11.691696, 13.165612, 14.601514, 13.070063, 
        14.315901, 15.534962, 16.747588, 17.980823, 12.852588, 13.891063, 
        14.855712, 15.636925, 16.977023, 9.6097598, 10.681909, 11.585263, 
        11.759212, 14.514652 };

        goodSymSab2 = { 1.6500137, 2.6179354, 3.1364259, 2.6765351, 1.3865851, 
        8.7834631, 8.2277943, 7.5389184, 5.4388668, 2.5003644, 13.070063, 
        11.501047, 10.013420, 6.9099146, 3.0690751, 12.852588, 11.163102, 
        9.5703208, 6.4333398, 2.8782114, 9.6097598, 8.5898752, 7.4599038, 
        4.8224008, 2.4451488 };

        THEN( "output scattering law vectors are correct" ){
          coldh( tev, ncold, transWeight + continWeight, 
            alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
          for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
            REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
            REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "molecular translations are assumed to be free" ){
      free = true;
      goodSymSab1 = { 0.639563665, 0.734847980, 0.642397321, 0.302707858, 
      1.452858E-2, 0.643890258, 0.756447167, 0.747722527, 0.573559847, 
      0.133577570, 0.731407289, 0.868992262, 0.913435993, 0.895184617, 
      0.537398952, 0.953382930, 1.139377425, 1.236241575, 1.375153242, 
      1.362943305, 5.380648E-2, 6.450705E-2, 7.117163E-2, 8.461831E-2, 
      0.108930287 };
      goodSymSab2 = { 0.6395636651, 0.5890740537, 0.4129826275, 0.1255578293, 
      2.4609421E-3, 0.6438902587, 0.6068176750, 0.4812008042, 0.2376524961, 
      2.2963400E-2, 0.7314072890, 0.6972235366, 0.5880200309, 0.3709822926, 
      9.2306067E-2, 0.9533829308, 0.9141851653, 0.7958588829, 0.5699187102, 
      0.2340743029, 5.3806485E-2, 5.1753352E-2, 4.5810754E-2, 3.5056533E-2, 
      1.8690471E-2 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
        }
      } // THEN
    } // WHEN
  } // GIVEN


  GIVEN( "Ortho deuterium"){
    ncold = 3;

    WHEN( "molecular translations are assumed to not be free" ){
      free = false;
      goodSymSab1 = { 1.04353549, 2.04283175, 3.04085742, 4.05944260, 
      5.05816298, 8.84191503, 10.3099118, 11.7696367, 13.2531960, 14.6990222, 
      27.5936205, 30.1836804, 32.7366103, 35.2069241, 37.9531925, 83.3216594, 
      89.4448115, 95.2689283, 99.3272795, 109.644553, 9.60975989, 10.6819098, 
      11.5852639, 11.7592124, 14.5146525 };

      goodSymSab2 = { 1.04353549, 1.64568291, 1.96673598, 1.67310377, 
      0.86654914, 8.84191503, 8.28251716, 7.58909537, 5.47516649, 2.51711392, 
      27.5936205, 24.2327449, 21.0875700, 14.5605371, 6.49983990, 83.3216594, 
      71.7845333, 61.3390811, 41.1194875, 18.7962242, 9.60975989, 8.58987526, 
      7.45990382, 4.82240089, 2.44514882 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
        }
      } // THEN
    } // WHEN
    WHEN( "molecular translations are assumed to be free" ){
      free = true;
      goodSymSab1 = { 0.639563665, 0.734847980, 0.642397321, 0.302707858, 
      1.452858E-2, 0.643890258, 0.756447167, 0.747722527, 0.573559847, 
      0.133577570, 0.731407289, 0.868992262, 0.913435993, 0.895184617, 
      0.537398952, 0.953382930, 1.139377425, 1.236241575, 1.375153242, 
      1.362943305, 5.380648E-2, 6.450705E-2, 7.117163E-2, 8.461831E-2, 
      0.108930287 };
      goodSymSab2 = { 0.6395636651, 0.5890740537, 0.4129826275, 0.1255578293, 
      2.4609421E-3, 0.6438902587, 0.6068176750, 0.4812008042, 0.2376524961, 
      2.2963400E-2, 0.7314072890, 0.6972235366, 0.5880200309, 0.3709822926, 
      9.2306067E-2, 0.9533829308, 0.9141851653, 0.7958588829, 0.5699187102, 
      0.2340743029, 5.3806485E-2, 5.1753352E-2, 4.5810754E-2, 3.5056533E-2, 
      1.8690471E-2 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
        }
      } // THEN
    } // WHEN
  } // GIVEN

  ska = { 0.001, 0.002, 0.005, 0.008, 0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0, 8.0, 5.0, 2.0, 1.0, 0.8, 0.5, 0.2, 0.1, 0.08, 0.05, 0.02, 0.01, 0.008, 0.005, 0.002, 0.001 };

  dka = 0.001;
  GIVEN( "Para deuterium"){
    ncold = 4;
    WHEN( "molecular translations are assumed to not be free" ){
      free = false;
      goodSymSab1 = { 1.501838, 2.955473, 4.393971, 5.790311, 7.232537, 
      7.989183, 9.355920, 10.68840, 11.95064, 13.19231, 12.11719, 13.37025, 
      14.58463, 15.67633, 16.58969, 12.40180, 13.52899, 14.75769, 15.74051, 
      16.42520, 9.698859, 10.68897, 12.42041, 13.60366, 15.10939 };

      goodSymSab2 = { 1.501838, 2.357724, 2.810726, 2.418721, 1.247006, 
      7.989183, 7.483196, 6.858405, 4.990571, 2.282003, 12.11719, 10.68925, 
      9.369173, 6.566440, 2.889771, 12.40180, 10.80009, 9.497284, 6.633682, 
      2.899509, 9.698859, 8.510247, 8.006278, 5.762811, 2.694623 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-5) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-5) );
        }
      } // THEN
    } // WHEN
    WHEN( "molecular translations are assumed to be free" ){
      free = true;
      goodSymSab1 = { 0.91727408, 1.0521462, 0.92011933, 0.43721466, 0.02316939,
      0.57958513, 0.68071185, 0.67378866, 0.52086685, 0.12924799, 0.32304386, 
      0.38388029, 0.40419732, 0.39909385, 0.24906591, 0.14814365, 0.17707755, 
      0.19234509, 0.21492490, 0.21699345, 0.05786461, 0.06936676, 0.07651663, 
      0.09089733, 0.11663430 };

      goodSymSab2 = { 0.91727408, 0.84528166, 0.59352393, 0.18101081, 
      0.00404511, 0.57958513, 0.54642996, 0.43411633, 0.21600343, 0.02215599, 
      0.32304386, 0.30805349, 0.26028263, 0.16546748, 0.04277909, 0.14814365, 
      0.14208293, 0.12383331, 0.08908330, 0.03727785, 0.05786461, 0.05565516, 
      0.04925662, 0.03766815, 0.02003144 };

      THEN( "output scattering law vectors are correct" ){
        coldh( tev, ncold, transWeight + continWeight, 
          alpha, beta, dka, ska, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
