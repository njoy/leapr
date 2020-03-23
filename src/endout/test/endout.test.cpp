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
  std::vector<double> alphas, betas, sab, temps, secondaryScatterVecThing (10,0.0), dwpix, dwp1;
  alphas = {1.1, 2.2, 3.3, 4.5, 5.8};
  betas = {0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7};
  sab = { 0.8812416231, 0.8226588624, 0.0371931173, 0.0370414517, 0.0453434977, 
          0.0449511778, 0.0177215328, 0.4944382555, 0.4922000032, 0.0833173024, 
          0.0740790498, 0.0698111325, 0.0693835223, 0.0339337297, 0.3223316052, 
          0.3256420103, 0.1180584841, 0.1055070878, 0.0815560280, 0.0810404205, 
          0.0458832427, 0.2182146674, 0.2226695868, 0.1282393202, 0.1182525679, 
          0.0860463732, 0.0851264533, 0.0537751752, 0.1498324933, 0.1534650775, 
          0.1195514020, 0.1135384304, 0.0852010403, 0.0844404310, 0.0572611451};
  temps = { 296.0 };
  dwpix = { 0.27366867553080776 };
  dwp1  = { 0.0 };






  std::cout.precision(15);
//std::cout << sab[0+4*betas.size()] << std::endl;
  double awr = 0.99917, spr = 20.449, aws = 15.85316, sps = 3.8883;

  int numSecondaryScatterers = 1, secondaryScatterType = 1;

  //endout(sab,awr,aws,spr,sps,temps,numSecondaryScatterers,secondaryScatterType,secondaryScatterVecThing,alphas,betas,dwpix,dwp1);

} // TEST CASE
