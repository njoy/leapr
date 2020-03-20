#include "catch.hpp"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "generalTools/testing.h"

std::string validSEND() {
  return
    "                                                                    27 7  0     \n";
}

std::string check1() {
  return
" 1.270000+2 8.934780+0          0          1          3          0  27 7  4      \n"
" 0.000000+0 0.000000+0          0          0         12          1  27 7  4      \n"
" 6.150000+0 6.000000-1 8.934780+0 1.518000-2 0.000000+0 1.000000+0  27 7  4      \n"
" 0.000000+0 3.748100+0 1.585800+1 0.000000+0 0.000000+0 1.000000+0  27 7  4      \n"
" 0.000000+0 0.000000+0          0          0          1          7  27 7  4      \n"
"          7          4                                              27 7  4      \n"
" 2.960000+2-6.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 1.066516-2 2.000000-1 2.082273-2 3.000000-1 2.930870-2  27 7  4      \n"
" 4.000000+2-6.000000-1          4          0          3          0  27 7  4      \n"
" 1.814988-2 3.526223-2 4.951282-2                                   27 7  4      \n"
" 1.200000+3-6.000000-1          4          0          3          0  27 7  4      \n"
" 1.318445-1 2.432061-1 3.293787-1                                   27 7  4      \n"
" 2.960000+2-4.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 1.418638-2 2.000000-1 2.592329-2 3.000000-1 3.452343-2  27 7  4      \n"
" 4.000000+2-4.000000-1          4          0          3          0  27 7  4      \n"
" 2.508672-2 4.497266-2 5.954545-2                                   27 7  4      \n"
" 1.200000+3-4.000000-1          4          0          3          0  27 7  4      \n"
" 1.634772-1 2.819640-1 3.704010-1                                   27 7  4      \n"
" 2.960000+2-2.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 2.137508-2 2.000000-1 3.892861-2 3.000000-1 5.227166-2  27 7  4      \n"
" 4.000000+2-2.000000-1          4          0          3          0  27 7  4      \n"
" 3.788786-2 6.838178-2 9.144674-2                                   27 7  4      \n"
" 1.200000+3-2.000000-1          4          0          3          0  27 7  4      \n"
" 2.700463-1 4.629903-1 5.986514-1                                   27 7  4      \n"
" 2.960000+2 0.000000+0          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 2.664996-2 2.000000-1 5.008647-2 3.000000-1 7.026869-2  27 7  4      \n"
" 4.000000+2 0.000000+0          4          0          3          0  27 7  4      \n"
" 4.833979-2 9.022656-2 1.254638-1                                   27 7  4      \n"
" 1.200000+3 0.000000+0          4          0          3          0  27 7  4      \n"
" 3.875195-1 6.638266-1 8.498093-1                                   27 7  4      \n"
" 2.960000+2 2.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 6.709659-3 2.000000-1 1.298148-2 3.000000-1 1.879842-2  27 7  4      \n"
" 4.000000+2 2.000000-1          4          0          3          0  27 7  4      \n"
" 1.258480-2 2.454353-2 3.554293-2                                   27 7  4      \n"
" 1.200000+3 2.000000-1          4          0          3          0  27 7  4      \n"
" 1.141221-1 2.132539-1 2.918076-1                                   27 7  4      \n"
" 2.960000+2 4.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 3.900015-3 2.000000-1 7.662907-3 3.000000-1 1.131221-2  27 7  4      \n"
" 4.000000+2 4.000000-1          4          0          3          0  27 7  4      \n"
" 7.544915-3 1.507582-2 2.236445-2                                   27 7  4      \n"
" 1.200000+3 4.000000-1          4          0          3          0  27 7  4      \n"
" 7.529282-2 1.467675-1 2.066500-1                                   27 7  4      \n"
" 2.960000+2 6.000000-1          2          0          1          3  27 7  4      \n"
"          3          4                                              27 7  4      \n"
" 1.000000-1 4.036312-3 2.000000-1 8.005619-3 3.000000-1 1.182832-2  27 7  4      \n"
" 4.000000+2 6.000000-1          4          0          3          0  27 7  4      \n"
" 8.044133-3 1.608329-2 2.370712-2                                   27 7  4      \n"
" 1.200000+3 6.000000-1          4          0          3          0  27 7  4      \n"
" 8.297012-2 1.572996-1 2.183293-1                                   27 7  4      \n"
" 0.000000+0 0.000000+0          0          0          1          3  27 7  4      \n"
"          3          2                                              27 7  4      \n"
" 2.960000+2 5.966722+2 4.000000+2 6.441872+2 1.200000+3 1.292336+3  27 7  4      \n"
" 0.000000+0 0.000000+0          0          0          1          3  27 7  4      \n"
"          3          2                                              27 7  4      \n"
" 2.960000+2 4.279227+2 4.000000+2 5.028527+2 1.200000+3 1.236600+3  27 7  4      \n";
}




TEST_CASE( "Preparing full ENDF output for S(a,b) --> [7,4]" ){
  GIVEN( "Be in BeO example - multiple temperatures and 1 secondary scatterer" ){
    using namespace njoy::ENDFtk;
    using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;

    std::vector<double> alphas         { 0.1, 0.2, 0.3 }, 
                        betas          { 0.0, 0.2 , 0.4, 0.6 },
                        temps          { 296.0, 400.0, 1200.0 },
                        effectiveTempsPrincipal { 596.6722, 644.18716, 1292.3357},
                        effectiveTempsSecondary { 427.9226, 502.85271, 1236.5996},

             //  S(a0,b0,T),  S(a0,b1,T), S(a0,b2,T),   S(a0,b3,T)
             //  S(a1,b0,T),  S(a1,b1,T), S(a1,b2,T),   S(a1,b3,T)
             //  S(a2,b0,T),  S(a2,b1,T), S(a2,b2,T),   S(a2,b3,T)
    sab_temp_1 {3.803356e-2, 1.186118e-2, 5.35523e-3,  5.494297e-3,
                7.283326e-2, 2.289232e-2, 1.153210e-2, 1.067055e-2,
                1.046095e-1, 3.313806e-2, 1.680395e-2, 1.554271e-2 }, 
    sab_temp_2 {6.888776e-2, 2.157749e-2, 1.085073e-2, 1.007578e-2,    
                1.308459e-1, 4.144943e-2, 2.104045e-2, 1.952161e-2,
                1.864122e-1, 5.972006e-2, 3.059757e-2, 2.836719e-2 },
    sab_temp_3 {5.744355e-1, 1.848110e-1, 9.506814e-2, 8.835550e-2, 
                1.011715e+0, 3.376009e-1, 1.787335e-1, 1.656053e-1,
                1.337622e+0, 4.627526e-1, 2.518537e-1, 2.328031e-1 };
    std::vector<std::vector<double>> fullSAB {sab_temp_1, sab_temp_2, sab_temp_3};

    int isym = 0, ilog = 0; // also known as lln

    unsigned int natoms_principal = 1;
    unsigned int natoms_secondary = 1;
    using std::move;
    double za  = 127.0, 
           awr_principal = 8.934780e+0, awr_secondary = 15.858, 
           xs_principal  = 6.153875e0,  xs_secondary  = 3.7481;

    int lasym = 0, lat = 1;


    //-------------- Create Scattering Law Constants Object -------------------
    int numSecondaryScatterers = 1;
    double epsilon = 1.976285e2, emax = 5.000001e0;
    std::vector<double>  
      xsVec  {move(xs_principal),  move(xs_secondary)  },
      awrVec {move(awr_principal), move(awr_secondary) };
    std::vector<unsigned int>  
      natomsVec{ move(natoms_principal), move(natoms_secondary) },
      secondaryScattererTypes { 0 }; 
      // Secondary Scatter Options: 0  = SCT, 1 = free, 2 = S(a,b)

    ScatteringLawConstants constants( ilog, numSecondaryScatterers, epsilon, 
      emax, std::move(xsVec), std::move(awrVec), std::move(natomsVec), 
      std::move(secondaryScattererTypes) );

    //-------------------------------------------------------------------------


    auto chunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                             effectiveTempsPrincipal, effectiveTempsSecondary, 
                             lasym, lat, isym, ilog, constants);
    std::string buffer;
    auto output = std::back_inserter( buffer );
    chunk.print( output, 27, 7 );
    //std::cout << buffer << std::endl;

    std::string sectionString = check1() + validSEND();

  } // GIVEN
} // TEST CASE


