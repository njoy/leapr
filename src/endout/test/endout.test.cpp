#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "endout/endout.h"
#include "generalTools/testing.h"
#include "coher/coher.h"


TEST_CASE( "finalize the debye-waller coefficient output" ){
  GIVEN( "No secondary scatterers considered" ){
  } // GIVEN
  GIVEN( "Secondary scatterers considered" ){
    WHEN( "SCT approximation (BeO test 23)" ){
        std::cout.precision(15);

        int numSecondaryScatterers = 1, secondaryScatterType = 0;
        double awr = 8.9347799999999999, aws = 15.858000000000001;
        std::vector<double> 

        temps { 296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0 },
        correctDWPIX {2.18024395, 2.72608260, 3.28130480, 3.85303452, 4.43480683, 
              5.02310000, 6.21177434, 7.41006538},
        correctDWP1 {2.16475875, 2.59570940, 3.04447306, 3.51490437, 3.99970405, 
              4.49435669, 5.50273307, 6.5268008},
        dwpix { 0.88189718, 1.49011626, 2.24201096, 3.15918678, 4.24222697, 
                5.49139855, 8.48861475, 12.1513474 },
        dwp1  { 0.49335305, 0.79941568, 1.17203004, 1.62375823, 2.15567111, 
                2.76830621, 4.23677157, 6.03029012 };

       scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awr, aws );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );


    } // WHEN
    WHEN( "Secondary scatterers use free gas approximation" ){
      THEN( "DWP1 is not changed and DWPIX is scaled by awr" ){
        int numSecondaryScatterers = 1, secondaryScatterType = 1;
        std::vector<double> 
        dwp1 (5,0.0),
        dwpix { 0.27366867, 0.27913484, 0.43494593, 0.86192089, 1.44790731 },
        temps { 296.0, 300.0, 400.0, 600.0, 800.0 },
        correctDWPIX{10.73794694, 10.80639075, 12.62883113, 16.68414797, 21.02028737 },
        correctDWP1 (5,0.0);
        
        double awr = 0.99917, aws = 15.85316; 
        scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awr, aws );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE






TEST_CASE( "endout" ){
  REQUIRE( true );
  std::vector<double> alphas, betas, sab, temps, secondaryScatterVecThing (10,0.0), dwpix, dwp1, tempf, tempf1;


  alphas = { 0.4, 4.0 };
  betas  = { 0.0, 2.0, 4.0, 8.0, 16.0 };
  sab    = { 0.202037135, 0.118385279, 0.019020055, 2.1675111E-3, 4.4194542E-6, 
             0.029581729, 0.062124313, 0.084467155, 0.0746431711, 0.0147378131 };
  temps  = { 293.6 };
  dwpix  = { 2.063703012 };
  dwp1   = { 1.893758364 };
  tempf  = { 486.4483635 };
  tempf1 = { 508.3412436 };
  double awr = 27.84423, spr = 2.021, aws = 15.862, sps = 7.4975;
  int numSecondaryScatterers = 1, secondaryScatterType = 0;
  int iel = 0;
  double translationalWeight = 0.0;
  std::cout.precision(15);
  std::vector<double> bragg(0);
  int numEdges = 0;
  int za = 147;
  int ilog = 0, isym = 0, lat = 1, numPrincipalAtoms = 1, numSecondaryAtoms = 1;


  endout(sab,za,awr,aws,spr,sps,temps,numSecondaryScatterers,secondaryScatterType,secondaryScatterVecThing,alphas,betas,dwpix,dwp1,iel,translationalWeight,bragg,numEdges,tempf,tempf1,ilog,isym,lat,numPrincipalAtoms,numSecondaryAtoms);

} // TEST CASE
