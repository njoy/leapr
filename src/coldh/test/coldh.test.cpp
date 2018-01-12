#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "discre/discre_util/bfill.h"
#include "discre/discre_util/exts.h"
#include "coldh/coldh.h"


void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-5 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-5 );
}

void equal_vec_mega_vec( std::vector<std::vector<std::vector<double>>> a, 
  std::vector<double> b ){
  REQUIRE( a.size()*a[0].size()*a[0][0].size() == b.size() );
  int i = 0;
  for ( auto a1 : a ){
    for ( auto a2 : a1 ){
      for ( auto a3 : a2 ){
        equal( a3, b[i] );
        i += 1;
      }
    }
  }
}



auto populateSymSab( const std::vector<double>& alpha, const std::vector<double>& beta, bool is_normal ){
  std::vector<std::vector<std::vector<double>>> sym_sab(alpha.size(),
    std::vector<std::vector<double>>(beta.size(),std::vector<double>(1,0.0)));
  int i = 1;
  for ( auto a = 0; a < alpha.size(); ++a ){
    for ( auto b = 0; b < beta.size(); ++b ){
      if ( is_normal ){ sym_sab[a][b][0] = i; }
      else {sym_sab[a][b][0] = 0; }
      i += 1;
    }
  }
  return sym_sab;
}



TEST_CASE( "coldh" ){

  std::vector<double> tempf {193093.99765}, alpha {0.1, 0.2, 0.4, 0.8, 1.6},
    beta {0.10, 0.15, 0.30, 0.60, 1.2}, ska { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13 },
    correctSymSab (25), correctSymSab2 (25);

  double temp = 200.0, tev = 1.723477e-2, tbeta = 2.0, 
    trans_weight = 0.3, scaling = 1.0, dka = 0.2;

  int itemp = 0, ncold = 1, nbeta = 5, lat = 3;

  bool free = false;

  auto sym_sab = populateSymSab( alpha, beta, true );
  auto sym_sab_2 = populateSymSab( alpha, beta, false );

  GIVEN( "molecular translations are assumed to not be free" ){

    WHEN( "step size is moderate" ){
      correctSymSab = { 1.7113874, 3.3919859, 5.0714083, 6.7348452, 8.4218850, 
        7.6149802, 8.8664244, 10.115396, 11.306691, 12.613772, 13.088162, 
        14.251586, 15.410194, 16.408195, 17.759561, 16.788298, 17.794618, 
        18.792099, 19.406237, 20.903111, 17.308492, 18.066429, 18.809744, 
        18.777265, 20.590641 };
      correctSymSab2 = { 1.7113874, 2.9207714, 3.7584709, 3.7033133, 2.5437978,
        7.6149802, 7.6447762, 7.5170871, 6.2697456, 3.8402507, 13.088162, 
        12.312841, 11.499779, 9.2237236, 5.4789194, 16.788298, 15.438241, 
        14.144248, 11.221266, 6.6252107, 17.308492, 15.813446, 14.417152, 
        11.526322, 6.8941688 };

      THEN( "output scattering law vectors are correct" ){
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling, 
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
      } // THEN
    } // WHEN


    WHEN( "step size is small" ){
      dka = 0.0001;

      correctSymSab = { 1.3287711, 2.6267533, 3.9235595, 5.2043802, 6.5088037, 
        7.6149802, 8.8664244, 10.115396, 11.306691, 12.613772, 13.088162, 
        14.251586, 15.410194, 16.408195, 17.759561, 16.788298, 17.794618, 
        18.792099, 19.406237, 20.903111, 17.308492, 18.066429, 18.809744, 
        18.777265, 20.590641 };

      correctSymSab2 = { 1.3287711, 2.2621297, 2.9081236, 2.8633763, 1.9675888,
        7.6149802, 7.6447762, 7.5170871, 6.2697456, 3.8402507, 13.088162, 
        12.312841, 11.499779, 9.2237236, 5.4789194, 16.788298, 15.438241, 
        14.144248, 11.221266, 6.6252107, 17.308492, 15.813446, 14.417152, 
        11.526322, 6.8941688 };

      THEN( "output scattering law vectors are correct" ){
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
      } // THEN

      sym_sab = populateSymSab( alpha, beta, true );
      sym_sab_2 = populateSymSab( alpha, beta, false );
 
      ska = { 0.1, 0.2, 0.3, 0.5, 0.8, 1.13 };

      THEN( "output scattering law vectors are correct" ){
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
      } // THEN
    } // WHEN


    WHEN( "step size is large" ){
      dka = 10.0;
      ska = { 100.001, 5.002, 0.003 };

      correctSymSab = { 4.712662, 9.394536, 14.07523, 18.73994, 23.42826, 
        26.36307, 30.73920, 35.11286, 39.42884, 43.86060, 43.00731, 46.89065, 
        50.76919, 54.48711, 58.55840, 50.71420, 53.84089, 56.95874, 59.69325, 
        63.31049, 45.48792, 47.58773, 49.67292, 50.98232, 54.137579 };

      correctSymSab2 = { 4.712662, 8.087215, 10.42866, 10.29185, 7.063631, 
        26.36307, 26.47085, 26.03566, 21.70350, 13.25161, 43.00731, 40.40555, 
        37.69436, 30.12187, 17.76729, 50.71420, 46.46355, 42.41879, 33.33124, 
        19.39806, 45.48792, 41.22266, 37.28116, 29.20083, 16.998312 };

      THEN( "output scattering law vectors are correct" ){
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "molecular translations are assumed to be free" ){
    ska = { 1.1, 2.2, 3.3, 4.4 };
    WHEN( "ortho hydrogen is chosen" ){
      tempf[0] = 200000; ncold = 1; dka = 2.3;
      free = true;

      correctSymSab = { 0.6906337, 0.7707476, 0.7716237, 0.6687087, 0.2834944, 
        0.4527882, 0.5086638, 0.5282524, 0.5302170, 0.4020220, 0.2736184, 
        0.3085437, 0.3264013, 0.3525404, 0.3571662, 0.1412989, 0.1596664, 
        0.1704781, 0.1909689, 0.2234552, 5.380649E-2, 6.086823E-2, 6.529124E-2, 
        7.447846E-2, 9.362654E-2 };
    
      correctSymSab2 = { 0.6906337, 0.6648041, 0.5741337, 0.3704882, 
        8.711346E-2,  0.4527882, 0.4388417, 0.3931755, 0.2936692, 0.1230125, 
        0.2736184, 0.2660678, 0.2427135, 0.1949027, 0.1090072, 0.1412989, 
        0.1376181, 0.1266453, 0.1053845, 6.801453E-2, 5.380649E-2, 5.244522E-2, 
        4.847128E-2, 4.104705E-2, 2.843402E-2 };

      THEN( "output scattering law vectors are correct" ){
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
      } // THEN
    } // WHEN

    WHEN( "para hydrogen is chosen" ){
      tempf[0] = 200000; ncold = 2; dka = 2.3;

      correctSymSab = { 3.761848E-2, 4.619964E-2, 5.011692E-2, 5.386082E-2, 
        4.389428E-2, 3.985561E-2, 4.787329E-2, 5.246895E-2, 6.015803E-2, 
        6.454883E-2, 4.580030E-2, 5.319663E-2, 5.767247E-2, 6.629526E-2, 
        7.942194E-2, 4.741764E-2, 5.420193E-2, 5.841111E-2, 6.698852E-2, 
        8.363218E-2, 3.549853E-2, 4.035129E-2, 4.343504E-2, 4.996043E-2, 
        6.417559E-2 };

      correctSymSab2 = { 3.761848E-2, 3.572692E-2, 2.999425E-2, 1.959557E-2, 
        8.295755E-3, 3.985561E-2, 3.826508E-2, 3.363161E-2, 2.537094E-2, 
        1.394832E-2, 4.580030E-2, 4.435546E-2, 4.013064E-2, 3.232050E-2, 
        1.984000E-2, 4.741764E-2, 4.610451E-2, 4.226779E-2, 3.511375E-2, 
        2.316111E-2, 3.549853E-2, 3.457250E-2, 3.188555E-2, 2.692709E-2, 
        1.866112E-2 };
    
      free = true;
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
 
    } // WHEN

    WHEN( "ortho deuterium is chosen" ){
      tempf[0] = 200000; ncold = 3; dka = 2.3;
      free = true;

      correctSymSab = { 1.380304, 1.543119, 1.546022, 1.340770, 0.5632300, 
        0.9730467, 1.095199, 1.138230, 1.142990, 0.8615740, 0.6178966, 
        0.6978264, 0.7386705, 0.7982229, 0.8065778, 0.3149662, 0.3563194, 
        0.3806761, 0.4268510, 0.4999492, 0.1121509, 0.1269906, 0.1363143, 
        0.1557608, 0.1967059 };

      correctSymSab2 = { 1.380304, 1.328098, 1.145199, 0.7357275, 0.1696729, 
        0.9730467, 0.9426315, 0.8431951, 0.6272528, 0.2594906, 0.6178966, 
        0.6006209, 0.5472136, 0.4380633, 0.2429255, 0.3149662, 0.3066842, 
        0.2820068, 0.2342521, 0.1505688, 0.1121509, 0.1092965, 0.1009743, 
        8.546640E-2, 5.922167E-2 };
      
    
      coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
 
    } // WHEN

    WHEN( "para deuterium is chosen" ){
      tempf[0] = 200000; ncold = 4; dka = 2.3;
      free = true;

      correctSymSab = { 1.288188, 1.439972, 1.442868, 1.252397, 0.5294604, 
        0.9198252, 1.035278, 1.076075, 1.081155, 0.8172565, 0.5941633, 
        0.6710235, 0.7103489, 0.7678541, 0.7769734, 0.3101028, 0.3508164, 
        0.3748090, 0.4203312, 0.4926102, 0.1138499, 0.1289094, 0.1383714, 
        0.1581061, 0.1996595 };

      correctSymSab2 = { 1.288188, 1.239549, 1.069146, 0.6875359, 0.1594064, 
        0.9198252, 0.8911030, 0.7972289, 0.5934155, 0.2461694, 0.5941633, 
        0.5775611, 0.5262499, 0.4214230, 0.2340327, 0.3101028, 0.3019515, 
        0.2776671, 0.2306855, 0.1483746, 0.1138499, 0.1109532, 0.1025076, 
        8.676988E-2, 6.013533E-2 };
      
    
      coldh( itemp, temp, tev, ncold, trans_weight, tbeta, tempf, scaling,
          alpha, beta, dka, ska, nbeta, lat, free, sym_sab, sym_sab_2 );
        equal_vec_mega_vec( sym_sab, correctSymSab );
        equal_vec_mega_vec( sym_sab_2, correctSymSab2 );
 
    } // WHEN
  } // GIVEN
} // TEST CASE
