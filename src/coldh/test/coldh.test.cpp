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

  double tbart = 2.1908278997651038; 
  lat = 1;

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
        coldh( tev, ncold, trans_weight, tbeta, scaling, 
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2, tbart );
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
        coldh( tev, ncold, trans_weight, tbeta, scaling, 
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2, tbart );
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
        coldh( tev, ncold, trans_weight, tbeta, scaling, 
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2, tbart );
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
        coldh( tev, ncold, trans_weight, tbeta, scaling, 
          alpha, beta, dka, ska, lat, free, sym_sab_1, sym_sab_2, tbart );
        for ( size_t i = 0; i < sym_sab_1.size(); ++i ){
          REQUIRE( goodSymSab1[i] == Approx(sym_sab_1[i]).epsilon(1e-6) );
          REQUIRE( goodSymSab2[i] == Approx(sym_sab_2[i]).epsilon(1e-6) );
        }
        REQUIRE(ranges::equal(sym_sab_1,goodSymSab1,equal));
        REQUIRE(ranges::equal(sym_sab_2,goodSymSab2,equal));
      } // THEN
    } // WHEN





  } // GIVEN
} // TEST CASE
