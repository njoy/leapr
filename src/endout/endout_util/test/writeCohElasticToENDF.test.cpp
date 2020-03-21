#include "catch.hpp"
#include "coher/coher.h"
#include "generalTools/testing.h"
#include "endout/endout_util/writeCohElasticToENDF.h"
#include "endout/endout_util/test/correctCoherentOutput.h"

using namespace njoy::ENDFtk;
using CoherentElastic = section::Type<7,2>::CoherentElastic;


TEST_CASE( "finalizing coherent elastic scattering data for ENDF" ){

  GIVEN( "one temperatures" ){
    WHEN( "no secondary scatterer present (aluminum)" ){
      int numSecondaryScatterers = 0, secondaryScatterType = 0;
      std::vector<double> bragg (60000,0.0);
      double maxEnergy = 5.0;
      int iel = 4, npr = 1;
      auto out = coher(iel,npr,bragg,maxEnergy);
      int numEdges = int(0.5*std::get<1>(out));

      double tol = 9e-8;
      std::vector<double> 
        temps { 20.0 },
        dwpix { 1.68005231 },
        dwp1  { 0.0 };

      auto myCohEl = writeCohElasticToENDF( bragg, dwpix, dwp1, 
           numSecondaryScatterers, secondaryScatterType, numEdges, tol, temps );

      THEN( "output CoherentElastic ENDF result is correct" ){
          auto begin = aluminum.begin();
          auto end = aluminum.end();
          long lineNumber = 1;
   
          CoherentElastic goodCohEl( begin, end, lineNumber, 101, 7, 2 );
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
      } // THEN
    } // WHEN
  } // GIVEN



  GIVEN( "three temperatures" ){
    WHEN( "secondary scatterer present (Beryllium Oxide)" ){
      int numSecondaryScatterers = 1, secondaryScatterType = 0;
      std::vector<double> bragg (10000,0.0);
      double maxEnergy = 5.0;
      int iel = 3, npr = 1;
      auto out = coher(iel,npr,bragg,maxEnergy);
      int numEdges = int(0.5*std::get<1>(out));

      double tol = 9e-8;
      std::vector<double> 
        temps { 296.0, 400.0, 1200.0 },
        dwpix { 2.1802439, 2.72608260, 7.41006538 },
        dwp1  { 2.1647587, 2.59570940, 6.52680080 };


      auto myCohEl = writeCohElasticToENDF( bragg, dwpix, dwp1, 
           numSecondaryScatterers, secondaryScatterType, numEdges, tol, temps );

      THEN( "output CoherentElastic ENDF result is correct" ){
          auto begin = berylliumOxide.begin();
          auto end = berylliumOxide.end();
          long lineNumber = 1;
   
          CoherentElastic goodCohEl( begin, end, lineNumber, 27, 7, 2 );
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

          REQUIRE( ranges::equal(temps,goodCohEl.temperatures(),equal) );
          REQUIRE( ranges::equal(temps, myCohEl.temperatures(),equal) );

          REQUIRE( ranges::equal(myCohEl.energies(),goodCohEl.energies(),equal) );

          for ( size_t itemp = 0; itemp < temps.size(); ++itemp ){
            auto correctXSVals = goodCohEl.thermalScatteringValues()[itemp];
            auto outputXSVals  =  myCohEl.thermalScatteringValues()[itemp];
            REQUIRE( ranges::equal(correctXSVals,outputXSVals,equal) );
          }
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "four temperatures" ){
    WHEN( "no secondary scatterer present (graphite)" ){
      int numSecondaryScatterers = 0, secondaryScatterType = 0;
      std::vector<double> bragg (10000,0.0);
      double maxEnergy = 1.0;
      int iel = 1, npr = 1;
      auto out = coher(iel,npr,bragg,maxEnergy);
      int numEdges = int(0.5*std::get<1>(out));

      double tol = 9e-8;
      std::vector<double> 
        temps { 296.0, 400.0, 600.0, 1200.0, 1400.0 },
        dwpix { 2.86028437, 3.63462793, 5.18435853,  10.00388006, 11.63097678 },
        dwp1  { 0.0, 0.0, 0.0, 0.0, 0.0 };


      auto myCohEl = writeCohElasticToENDF( bragg, dwpix, dwp1, 
           numSecondaryScatterers, secondaryScatterType, numEdges, tol, temps );

      THEN( "output CoherentElastic ENDF result is correct" ){
          auto begin = graphite.begin();
          auto end = graphite.end();
          long lineNumber = 1;
   
          CoherentElastic goodCohEl( begin, end, lineNumber, 30, 7, 2 );
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

          REQUIRE( ranges::equal(temps,goodCohEl.temperatures(),equal) );
          REQUIRE( ranges::equal(temps,  myCohEl.temperatures(),equal) );

          REQUIRE( ranges::equal(myCohEl.energies(),goodCohEl.energies(),equal) );

          for ( size_t itemp = 0; itemp < temps.size(); ++itemp ){
            auto correctXSVals = goodCohEl.thermalScatteringValues()[itemp];
            auto outputXSVals  =  myCohEl.thermalScatteringValues()[itemp];
            REQUIRE( ranges::equal(correctXSVals,outputXSVals,equal) );
          }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE

