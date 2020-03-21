#include "catch.hpp"
#include "endout/endout_util/writeInelasticToENDF.h"
#include "generalTools/testing.h"
#include "endout/endout_util/test/correctInelasticOutput.h"



TEST_CASE( "Preparing full ENDF output for S(a,b) --> [7,4]" ){
  GIVEN( "Be in BeO example - multiple temperatures and 1 secondary scatterer" ){
    using namespace njoy::ENDFtk;
    using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
    using Inelastic = section::Type< 7, 4 >;

    std::vector<double> alphas         { 0.1, 0.2, 0.3 }, 
                        betas          { 0.0, 0.2 , 0.4, 0.6 },
                        temps          { 296.0, 400.0, 1200.0 },
                        effectiveTempsPrincipal { 596.6722, 644.18716, 1292.3357},
                        effectiveTempsSecondary { 427.9226, 502.85271, 1236.5996},

             //  S(a0,b0,T),  S(a0,b1,T), S(a0,b2,T),   S(a0,b3,T)
             //  S(a1,b0,T),  S(a1,b1,T), S(a1,b2,T),   S(a1,b3,T)
             //  S(a2,b0,T),  S(a2,b1,T), S(a2,b2,T),   S(a2,b3,T)
    /*sab_temp_1 {3.803356e-2, 1.186118e-2, 5.35523e-3,  5.494297e-3,
                7.283326e-2, 2.289232e-2, 1.153210e-2, 1.067055e-2,
                1.046095e-1, 3.313806e-2, 1.680395e-2, 1.554271e-2 }, 
    sab_temp_2 {6.888776e-2, 2.157749e-2, 1.085073e-2, 1.007578e-2,    
                1.308459e-1, 4.144943e-2, 2.104045e-2, 1.952161e-2,
                1.864122e-1, 5.972006e-2, 3.059757e-2, 2.836719e-2 },
    sab_temp_3 {5.744355e-1, 1.848110e-1, 9.506814e-2, 8.835550e-2, 
                1.011715e+0, 3.376009e-1, 1.787335e-1, 1.656053e-1,
                1.337622e+0, 4.627526e-1, 2.518537e-1, 2.328031e-1 };
                */

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

    int isym = 0, ilog = 0; // also known as lln

    unsigned int natoms_principal = 1;
    unsigned int natoms_secondary = 1;
    using std::move;
    double za  = 127.0, 
           awr_principal = 8.934780e+0, awr_secondary = 15.858, 
           xs_principal  = 6.153875,  xs_secondary  = 3.7481;

    int lasym = 0, lat = 1;


    //-------------- Create Scattering Law Constants Object -------------------
    int numSecondaryScatterers = 1;
    double epsilon = betas[betas.size()-1], emax = 0.0253*betas[betas.size()-1];//5.000001e0;
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


    auto chunk1 = writeInelasticToENDF(fullSAB, alphas, betas, temps, za, 
                              effectiveTempsPrincipal, effectiveTempsSecondary, 
                              lasym, lat, isym, ilog, constants);


    auto begin = ENDF_BeO_correct.begin();
    auto end = ENDF_BeO_correct.end();
    long lineNumber = 1;
    HeadRecord head( begin, end, lineNumber );

    Inelastic chunk2( head, begin, end, lineNumber, 27 );


    REQUIRE( chunk1.ZA() == Approx(chunk2.ZA()) );
    REQUIRE( chunk1.AWR() == Approx(chunk2.AWR()) );

    REQUIRE( chunk1.LAT() == chunk2.LAT() );
    REQUIRE( chunk1.temperatureOption() == chunk2.temperatureOption() );
    REQUIRE( chunk1.LASYM() == chunk2.LASYM() );
    REQUIRE( chunk1.symmetryOption() == chunk2.symmetryOption() );

    auto barray1 = chunk1.constants();
    auto barray2 = chunk2.constants();
    REQUIRE( barray1.LLN() == barray2.LLN() );
    REQUIRE( barray1.sabStorageType() == barray2.sabStorageType() );
    REQUIRE( barray1.NI() == barray2.NI() );
    REQUIRE( barray1.numberConstants() == barray2.numberConstants() );
    REQUIRE( barray1.NS() == barray2.NS() );
    REQUIRE( barray1.numberNonPrincipalScatterers() == 
             barray2.numberNonPrincipalScatterers() );

    
    REQUIRE( barray2.epsilon() == Approx( barray2.epsilon() ) );
    REQUIRE( barray1.upperEnergyLimit() == Approx( barray2.upperEnergyLimit() ) );
    REQUIRE( ranges::equal(barray1.totalFreeCrossSections(),
                           barray2.totalFreeCrossSections(), equal) );
    REQUIRE( ranges::equal(barray1.atomicWeightRatios(),
                           barray2.atomicWeightRatios(), equal) );
    REQUIRE( ranges::equal(barray1.numberAtoms(),
                           barray2.numberAtoms(), equal) );
    REQUIRE( ranges::equal(barray1.analyticalFunctionTypes(),
                           barray2.analyticalFunctionTypes(), equal) );
  
    using Tabulated = section::Type< 7, 4 >::Tabulated;
  
    auto table1 = std::get< Tabulated >( chunk1.scatteringLaw() );
    auto table2 = std::get< Tabulated >( chunk2.scatteringLaw() );
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
  
    auto temp1 = chunk1.principalEffectiveTemperature();
    auto temp2 = chunk2.principalEffectiveTemperature();
    REQUIRE( temp1.NT() == temp2.NT() );
    REQUIRE( temp1.numberTemperatures() == temp2.numberTemperatures() );
    REQUIRE( temp1.NR() == temp2.NR() );
    REQUIRE( ranges::equal(temp1.interpolants(), temp2.interpolants(), equal) );
    REQUIRE( ranges::equal(temp1.boundaries(), temp2.boundaries(), equal) );
    REQUIRE( ranges::equal(temp1.moderatorTemperatures(), temp2.moderatorTemperatures(), equal) );
    REQUIRE( ranges::equal(temp1.effectiveTemperatures(), temp2.effectiveTemperatures(), equal) );

  
    std::string buffer;
    auto output = std::back_inserter( buffer );
    chunk1.print( output, 27, 7 );
    //std::cout << buffer << std::endl;
    //std::cout << std::endl;
    

/*
  // one secondary scatterer => 1 secondary temperature (std::nullopt in this
  // case)
  REQUIRE( 1 == chunk.secondaryEffectiveTemperatures().size() );
  REQUIRE( std::nullopt == chunk.secondaryEffectiveTemperatures()[0] );

  REQUIRE( 21 == chunk.NC() );
 
 
*/

    /*
    REQUIRE( chunk1.() == chunk2.() );
    REQUIRE( chunk1.() == chunk2.() );


    REQUIRE( chunk1.() == Approx(chunk2.()) );
    REQUIRE( chunk1.() == Approx(chunk2.()) );
    REQUIRE( chunk1.() == Approx(chunk2.()) );
    REQUIRE( chunk1.() == Approx(chunk2.()) );
    REQUIRE( chunk1.() == Approx(chunk2.()) );
    */


  } // GIVEN
} // TEST CASE


