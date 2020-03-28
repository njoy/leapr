#include <string>

using namespace njoy::ENDFtk;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
using Elastic = section::Type<7,2>;
using Inelastic = section::Type<7,4>;




inline void checkFullCohElastic(CoherentElastic goodCohEl, CoherentElastic myCohEl, const std::vector<double>& temps){
  REQUIRE( myCohEl.elasticScatteringType() == 
           goodCohEl.elasticScatteringType() );
  REQUIRE( myCohEl.temperatureDependenceFlag() == 
           goodCohEl.temperatureDependenceFlag() );
  REQUIRE( myCohEl.temperatureDependenceFlag() == 
           goodCohEl.temperatureDependenceFlag() );

  REQUIRE( myCohEl.NT() == goodCohEl.NT() );
  REQUIRE( myCohEl.NP() == goodCohEl.NP() );
  REQUIRE( myCohEl.NR() == goodCohEl.NR() );
  REQUIRE( myCohEl.NC() == goodCohEl.NC() );

  REQUIRE( myCohEl.numberBraggEdges() == goodCohEl.numberBraggEdges() );
  REQUIRE( myCohEl.numberBraggEdges() == int(myCohEl.energies().size()) );

  REQUIRE( ranges::equal(myCohEl.LI(), goodCohEl.LI(), equal) );

  REQUIRE( ranges::equal(goodCohEl.boundaries(), 
                         myCohEl.boundaries(), equal) );
  REQUIRE( ranges::equal(goodCohEl.interpolants(),
                         myCohEl.interpolants(), equal) );

  REQUIRE( ranges::equal(temps, goodCohEl.temperatures(), equal) );
  REQUIRE( ranges::equal(temps,   myCohEl.temperatures(), equal) );

  REQUIRE( ranges::equal(myCohEl.energies(),goodCohEl.energies(),equal) );

  for ( size_t itemp = 0; itemp < temps.size(); ++itemp ){
    auto correctXSVals = goodCohEl.thermalScatteringValues()[itemp];
    auto outputXSVals  =   myCohEl.thermalScatteringValues()[itemp];
    REQUIRE( ranges::equal(correctXSVals,outputXSVals,equal) );
  }


}


template <typename IncoherentElastic >
inline void testIncoherentElasticOutput( IncoherentElastic mine, IncoherentElastic good ){
  REQUIRE( good.LTHR() == mine.LTHR() );
  REQUIRE( good.elasticScatteringType() == mine.elasticScatteringType() );
  REQUIRE( good.SB() == Approx( mine.SB() ) );
  REQUIRE( good.boundCrossSection() == Approx( mine.boundCrossSection() ) );
  REQUIRE( good.NP() == mine.NP() );
  REQUIRE( good.numberTemperatures() == mine.numberTemperatures() );
  REQUIRE( good.NR() == mine.NR() );
  REQUIRE( ranges::equal(good.interpolants(), mine.interpolants(), equal) );
  REQUIRE( ranges::equal(good.boundaries(), mine.boundaries(), equal) );
  REQUIRE( ranges::equal(good.temperatures(), mine.temperatures(), equal) );
  REQUIRE( ranges::equal(good.debyeWallerValues(), mine.debyeWallerValues(), equal) );
}



template <typename ENDFtk_Inelastic>
inline void checkFullInelastic(std::string correctString, ENDFtk_Inelastic testChunk,
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


