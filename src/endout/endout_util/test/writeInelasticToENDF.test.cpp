#include "catch.hpp"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "generalTools/testing.h"
#include "endout/endout_util/test/correctInelasticOutput.h"
#include "endout/endout_util/test/check_MF_Output.h"

using namespace njoy::ENDFtk;
using ScatteringLawConstants = section::Type<7,4>::ScatteringLawConstants;
using Inelastic = section::Type<7,4>;


TEST_CASE( "Preparing full ENDF output for S(a,b) --> [7,4]" ){
  GIVEN( "Secondary scatterer example (Be in BeO)" ){

    int isym = 0, ilog = 0; // also known as lln

    unsigned int nAtomsPrincipal = 1;
    unsigned int nAtomsSecondary = 1;
    using std::move;
    double za  = 127.0, 
           awr_principal = 8.934780e+0, awr_secondary = 15.858, 
           xs_principal  = 6.153875,  xs_secondary  = 3.7481;
    int lasym = 0, lat = 1;

    std::vector<double> alphas         { 0.1, 0.2, 0.3 }, 
                        betas          { 0.0, 0.2 , 0.4, 0.6 };

    //-------------- Create Scattering Law Constants Object -------------------
    int numSecondaryScatterers = 1;
    double epsilon = betas[betas.size()-1], emax = 0.0253*betas[betas.size()-1];

    std::vector<double>  
      xsVec  {move(xs_principal),  move(xs_secondary)  },
      awrVec {move(awr_principal), move(awr_secondary) };
    std::vector<unsigned int>  
      natomsVec{ move(nAtomsPrincipal), move(nAtomsSecondary) },
      secondaryScattererTypes { 0 }; // 0 = SCT, 1 = Free, 2 = S(a,b)

    ScatteringLawConstants constants( ilog, numSecondaryScatterers, epsilon, 
      emax, std::move(xsVec), std::move(awrVec), std::move(natomsVec), 
      std::move(secondaryScattererTypes) );


    WHEN( "one temperature is considered" ){
      std::vector<double> temps          { 800.0 },
                          effectiveTempsPrincipal { 935.4821895 },
                          effectiveTempsSecondary { 854.2996689 },

               //  S(a0,b0,T),  S(a0,b1,T), S(a0,b2,T),   S(a0,b3,T)
      sab_temp_1 { 0.2655523746, 0.0874977193, 0.0462063716, 0.0445819048, 
                   0.4863204544, 0.1641521631, 0.0883833504, 0.0850936470, 
                   0.6682223599, 0.2310164432, 0.1267533036, 0.1218130094 };
      std::vector<std::vector<double>> fullSAB {sab_temp_1};

      auto myChunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                     effectiveTempsPrincipal, effectiveTempsSecondary, lasym, 
                     lat, isym, ilog, constants);

      auto begin = ENDF_BeO_correct_1temp.begin();
      auto end   = ENDF_BeO_correct_1temp.end();
      long lineNumber = 1;
      HeadRecord head( begin, end, lineNumber );
      Inelastic trueChunk( head, begin, end, lineNumber, 27 );
      checkFullInelastic( trueChunk, myChunk, betas );

    } // WHEN
    WHEN( "three temperatures are considered" ){
      std::vector<double> temps                   { 296.0, 400.0, 1200.0 },
                          effectiveTempsPrincipal { 596.6722, 644.18716, 1292.3357},
                          effectiveTempsSecondary { 427.9226, 502.85271, 1236.5996},

               //  S(a0,b0,T),  S(a0,b1,T), S(a0,b2,T),   S(a0,b3,T)
               //  S(a1,b0,T),  S(a1,b1,T), S(a1,b2,T),   S(a1,b3,T)
               //  S(a2,b0,T),  S(a2,b1,T), S(a2,b2,T),   S(a2,b3,T)
      sab_temp_1 { 0.0380264333, 0.0130957910, 0.0072363504, 0.0073965617, 
                   0.0728195676, 0.0252751190, 0.0140594421, 0.0143649433, 
                   0.1045897671, 0.0365872656, 0.0204866250, 0.0209239351 },
      sab_temp_2 { 0.0688748433, 0.0232169324, 0.0125637553, 0.0125544215, 
                   0.1308213751, 0.0445986629, 0.0243620961, 0.0243238989, 
                   0.1863772584, 0.0642573487, 0.0354279136, 0.0353454291 },
      sab_temp_3 { 0.5743284250, 0.1893563057, 0.0998141316, 0.0950600117, 
                   1.0115272080, 0.3459026996, 0.1876553149, 0.1781711799, 
                   1.3373756330, 0.4741303254, 0.2644242176, 0.2504673594 };

      std::vector<std::vector<double>> fullSAB {sab_temp_1, sab_temp_2, sab_temp_3};

      AND_WHEN( "No log option chosen (return S(a,b))" ){
        ilog = 0;
        auto myChunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                       effectiveTempsPrincipal, effectiveTempsSecondary, lasym, 
                       lat, isym, ilog, constants);
        auto begin = ENDF_BeO_correct_3temps.begin();
        auto end   = ENDF_BeO_correct_3temps.end();
        long lineNumber = 1;
        HeadRecord head( begin, end, lineNumber );
        Inelastic trueChunk( head, begin, end, lineNumber, 27 );
        checkFullInelastic( trueChunk, myChunk, betas );
      } // AND WHEN
      AND_WHEN( "No log option chosen (return S(a,b))" ){
        ilog = 1;
        std::vector<double>  xsVec  {move(xs_principal),  move(xs_secondary) },
                             awrVec {move(awr_principal), move(awr_secondary)};
        std::vector<unsigned int>  
          natomsVec{ move(nAtomsPrincipal), move(nAtomsSecondary) },
          secondaryScattererTypes { 0 }; // 0 = SCT, 1 = Free, 2 = S(a,b)
        ScatteringLawConstants constants( ilog, numSecondaryScatterers, epsilon, 
          emax, std::move(xsVec), std::move(awrVec), std::move(natomsVec), 
          std::move(secondaryScattererTypes) );
        auto myChunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                       effectiveTempsPrincipal, effectiveTempsSecondary, lasym, 
                       lat, isym, ilog, constants);
        auto begin = ENDF_BeO_correct_3temps_log.begin();
        auto end   = ENDF_BeO_correct_3temps_log.end();
        long lineNumber = 1;
        HeadRecord head( begin, end, lineNumber );
        Inelastic trueChunk( head, begin, end, lineNumber, 27 );
        checkFullInelastic( trueChunk, myChunk, betas );
      } // AND WHEN
      AND_WHEN( "Symmetric scattering law requested" ){
          /*
           * DOES NOT WORK BECAUSE ENDFTK NOT EQUIPPED TO READ LEAPR FILES WHEN
           * WE HAVE ISYM IN THERE (NOT LASYM)
        ilog = 0;
        isym = 2;
        std::vector<double>  xsVec  {move(xs_principal),  move(xs_secondary) },
                             awrVec {move(awr_principal), move(awr_secondary)};
        std::vector<unsigned int>  
          natomsVec{ move(nAtomsPrincipal), move(nAtomsSecondary) },
          secondaryScattererTypes { 0 }; // 0 = SCT, 1 = Free, 2 = S(a,b)
        ScatteringLawConstants constants( ilog, numSecondaryScatterers, epsilon, 
          emax, std::move(xsVec), std::move(awrVec), std::move(natomsVec), 
          std::move(secondaryScattererTypes) );
        auto myChunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                       effectiveTempsPrincipal, effectiveTempsSecondary, lasym, 
                       lat, isym, ilog, constants);
        checkFullInelastic( ENDF_BeO_correct_3temps_nonsymmetric, myChunk, betas );
        */
      } // AND WHEN

    } // WHEN
  } // GIVEN
} // TEST CASE


