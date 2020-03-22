#include "catch.hpp"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "generalTools/testing.h"
#include "endout/endout_util/test/correctInelasticOutput.h"

using namespace njoy::ENDFtk;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using Inelastic = section::Type< 7, 4 >;



template <typename ENDFtk_Inelastic>
void checkFullInelastic(std::string correctString, ENDFtk_Inelastic testChunk,
        std::vector<double> betas ){

    auto begin = correctString.begin();
    auto end = correctString.end();
    long lineNumber = 1;
    HeadRecord head( begin, end, lineNumber );

    Inelastic trueChunk( head, begin, end, lineNumber, 27 );

    REQUIRE( testChunk.ZA()  == Approx(trueChunk.ZA())  );
    REQUIRE( testChunk.AWR() == Approx(trueChunk.AWR()) );

    REQUIRE( testChunk.NC()    == trueChunk.NC()    );
    REQUIRE( testChunk.LAT()   == trueChunk.LAT()   );
    REQUIRE( testChunk.LASYM() == trueChunk.LASYM() );
    REQUIRE( testChunk.temperatureOption() == trueChunk.temperatureOption() );
    REQUIRE( testChunk.symmetryOption()    == trueChunk.symmetryOption()    );

    auto constsMine = testChunk.constants();
    auto constsNjoy = trueChunk.constants();
    REQUIRE( constsMine.LLN() == constsNjoy.LLN() );
    REQUIRE( constsMine.NI()  == constsNjoy.NI()  );
    REQUIRE( constsMine.NS()  == constsNjoy.NS()  );
    REQUIRE( constsMine.sabStorageType()  == constsNjoy.sabStorageType()  );
    REQUIRE( constsMine.numberConstants() == constsNjoy.numberConstants() );
    REQUIRE( constsMine.numberNonPrincipalScatterers() == 
             constsNjoy.numberNonPrincipalScatterers() );
    
    REQUIRE( constsMine.epsilon()          == Approx( constsNjoy.epsilon() ) );
    REQUIRE( constsMine.upperEnergyLimit() == Approx( constsNjoy.upperEnergyLimit() ) );

    REQUIRE( ranges::equal(constsMine.totalFreeCrossSections(),
                           constsNjoy.totalFreeCrossSections(), equal) );
    REQUIRE( ranges::equal(constsMine.atomicWeightRatios(),
                           constsNjoy.atomicWeightRatios(), equal) );
    REQUIRE( ranges::equal(constsMine.numberAtoms(),
                           constsNjoy.numberAtoms(), equal) );
    REQUIRE( ranges::equal(constsMine.analyticalFunctionTypes(),
                           constsNjoy.analyticalFunctionTypes(), equal) );
  
    using Tabulated = section::Type< 7, 4 >::Tabulated;
  
    auto table1 = std::get< Tabulated >( testChunk.scatteringLaw() );
    auto table2 = std::get< Tabulated >( trueChunk.scatteringLaw() );
    REQUIRE( table1.NR() == table2.NR() );
    REQUIRE( table1.NB() == table2.NB() );
    REQUIRE( table1.numberBetas() == table2.numberBetas() );
    REQUIRE( ranges::equal(table1.boundaries(),table2.boundaries(),equal) );
    REQUIRE( ranges::equal(table1.interpolants(),table2.interpolants(),equal) );
  
    for (size_t b = 0; b < betas.size(); ++b){
      auto value1 = table1.betas()[b];
      auto value2 = table2.betas()[b];
      REQUIRE( value1.beta() == Approx( value2.beta() ) );
      REQUIRE( value1.LT() == value2.LT() );
      REQUIRE( value1.temperatureDependenceFlag() == value2.temperatureDependenceFlag() );
      REQUIRE( value1.NT() == value2.NT() );
      REQUIRE( value1.numberTemperatures() == value2.numberTemperatures() );
  
      REQUIRE( value1.NR() == value2.NR() );
      REQUIRE( value1.NA() == value2.NA() );
      REQUIRE( value1.numberAlphas() == value2.numberAlphas() );
      REQUIRE( ranges::equal(value1.boundaries(),value2.boundaries(),equal) );
      REQUIRE( ranges::equal(value1.interpolants(),value2.interpolants(),equal) );
  

      REQUIRE( ranges::equal(value1.temperatures(), value2.temperatures(), equal) );
      REQUIRE( ranges::equal(value1.alphas(), value2.alphas(), equal) );
      REQUIRE( ranges::equal(value1.LI(), value2.LI(), equal) );
      REQUIRE( ranges::equal(value1.temperatureInterpolants(), 
                             value2.temperatureInterpolants(), 
                             [](auto x, auto y){return x == y;} ) );
      REQUIRE( value1.thermalScatteringValues().size() == 
               value2.thermalScatteringValues().size() );
      for (size_t i = 0; i < value1.thermalScatteringValues().size(); ++i){
        REQUIRE( ranges::equal(value1.thermalScatteringValues()[i], 
                               value2.thermalScatteringValues()[i], equal) );
      }
    }
  
    auto tempMine = testChunk.principalEffectiveTemperature();
    auto tempNjoy = trueChunk.principalEffectiveTemperature();
    REQUIRE( tempMine.NT() == tempNjoy.NT() );
    REQUIRE( tempMine.NR() == tempNjoy.NR() );
    REQUIRE( tempMine.numberTemperatures() == tempNjoy.numberTemperatures() );
    REQUIRE( ranges::equal(tempMine.interpolants(), tempNjoy.interpolants(), equal) );
    REQUIRE( ranges::equal(tempMine.boundaries(), tempNjoy.boundaries(), equal) );
    REQUIRE( ranges::equal(tempMine.moderatorTemperatures(), 
                           tempNjoy.moderatorTemperatures(), equal) );
    REQUIRE( ranges::equal(tempMine.effectiveTemperatures(), 
                           tempNjoy.effectiveTemperatures(), equal) );

  
    std::string buffer;
    auto output = std::back_inserter( buffer );
    testChunk.print( output, 27, 7 );
    }



TEST_CASE( "Preparing full ENDF output for S(a,b) --> [7,4]" ){
  GIVEN( "Secondary scatterer example (Be in BeO)" ){

    int isym = 0, ilog = 0; // also known as lln

    unsigned int natoms_principal = 1;
    unsigned int natoms_secondary = 1;
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
      natomsVec{ move(natoms_principal), move(natoms_secondary) },
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
      checkFullInelastic( ENDF_BeO_correct_1temp, myChunk, betas );

    } // WHEN
    WHEN( "three temperatures are considered" ){

      std::vector<double> temps          { 296.0, 400.0, 1200.0 },
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

      auto myChunk = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                     effectiveTempsPrincipal, effectiveTempsSecondary, lasym, 
                     lat, isym, ilog, constants);
      checkFullInelastic( ENDF_BeO_correct_3temps, myChunk, betas );


    } // WHEN
  } // GIVEN
} // TEST CASE


