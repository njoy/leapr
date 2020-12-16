#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "endout/endout.h"
#include "generalTools/testing.h"
#include "coherentElastic/coherentElastic.h"
#include "endout/test/correctMF7.h"
#include "endout/endout_util/test/check_MF_Output.h"
#include "endout/test/braggVectors.h"

template <typename T> std::string type_name();
using namespace njoy::ENDFtk;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
using MF7 = njoy::ENDFtk::file::Type<7>;

/*

void checkFullMF7( MF7 myMF7, MF7 goodMF7, const std::vector<double>& betas ){

    njoy::ENDFtk::section::Type<7,2> myMT2 = myMF7.MT(2_c);
    njoy::ENDFtk::section::Type<7,4> myMT4 = myMF7.MT(4_c);

    REQUIRE( goodMF7.MF()     == myMF7.MF() ); 
    REQUIRE( goodMF7.hasMT(2) == myMF7.hasMT(2) ); 
    REQUIRE( goodMF7.hasMT(4) == myMF7.hasMT(4) ); 


    njoy::ENDFtk::section::Type<7,4> g4 = goodMF7.MT(4_c), m4 = myMF7.MT(4_c);
    REQUIRE( g4.ZA()   == m4.ZA() );
    checkFullInelastic(g4,m4,betas);


    njoy::ENDFtk::section::Type<7,2> g2 = goodMF7.MT(2_c), m2 = myMF7.MT(2_c);
    REQUIRE( g2.ZA()   == m2.ZA() );
    REQUIRE( g2.LTHR() == m2.LTHR() );
    REQUIRE( g2.NC()   == m2.NC() );


    REQUIRE( g2.elasticScatteringType() == m2.elasticScatteringType() );
    if (g2.elasticScatteringType() == 2){
      auto g2_law = std::get<IncoherentElastic>(g2.scatteringLaw());
      auto m2_law = std::get<IncoherentElastic>(m2.scatteringLaw());
      testIncoherentElasticOutput(g2_law,m2_law);
    }
    //else if ( ){
    //checkFullCohElastic(g2_law,m2_law,temps);
    //}
}
    


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
                2.76830621, 4.23677157, 6.03029012 },
        awrVec { awr, aws };

       scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awrVec );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );


    } // WHEN
    WHEN( "Secondary scatterers use free gas approximation" ){
      THEN( "DWP1 is not changed and DWPIX is scaled by awr" ){
        int numSecondaryScatterers = 1, secondaryScatterType = 1;
        double awr = 0.99917, aws = 15.85316; 
        std::vector<double> 
        dwp1 (5,0.0),
        dwpix { 0.27366867, 0.27913484, 0.43494593, 0.86192089, 1.44790731 },
        temps { 296.0, 300.0, 400.0, 600.0, 800.0 },
        correctDWPIX{10.73794694, 10.80639075, 12.62883113, 16.68414797, 21.02028737 },
        correctDWP1 (5,0.0),
        awrVec { awr, aws };
        
        scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awrVec );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE






TEST_CASE( "endout" ){
  std::cout.precision(15);
  std::vector<double> alphas, betas, temps, dwpix, dwp1, tempf, tempf1;
  std::vector<std::vector<double>> sab, principalScatterSAB;

  GIVEN( "H2O example from ENDFB-VIII.0 library" ){
    WHEN( "No incoherent scattering and Inelastic Scattering" ){
      alphas = { 1e-5, 1e-3, 2.0, 12.0 };
      betas  = { 0.0,  1e-5, 1.0, 9.0 };
      sab    = { { 618568.32494, 1524.6039481, 7.3117729E-7, 2.348787E-14, 
                   6176.3232863, 5934.9644697, 7.3128682E-5, 3.459404E-10, 
                   0.7657006321, 0.7657019792, 0.1022689961, 2.8459119E-3, 
                   2.4934839E-2, 2.4934968E-2, 3.3381878E-2, 2.6562263E-2 }, 
                 { 62683.075534, 29788.408139, 5.4697371E-6, 1.335162E-13, 
                   631.28974258, 629.99613584, 5.4808088E-4, 1.9584902E-9, 
                   0.4220947953, 0.4220958818, 0.3427480220, 8.7691637E-3, 
                   4.4382878E-2, 4.4382997E-2, 5.4603240E-2, 5.2083012E-2 }, 
                 { 50853.269587, 31557.012520, 7.8831344E-6, 1.825618E-13, 
                   514.57830261, 513.82421213, 7.8951831E-4, 2.6961786E-9, 
                   0.4020187562, 0.4020197225, 0.3648925394, 9.4606148E-3, 
                   4.8354693E-2, 4.8354811E-2, 5.8593715E-2, 5.6786087E-2 }, 
                 { 41980.988868, 30941.069562, 1.0404665E-5, 2.329200E-13, 
                   427.13820831, 426.79979962, 1.0433031E-3, 3.4550042E-9, 
                   0.3986884832, 0.3986893708, 0.3805113741, 1.0098559E-2, 
                   5.2357066E-2, 5.2357184E-2, 6.2662008E-2, 6.1317788E-2 }, 
                 { 51665.912938, 38098.033371, 1.5150875E-5, 3.465907E-13, 
                   522.87835724, 522.54084261, 1.5194978E-3, 5.1423444E-9, 
                   0.4311205135, 0.4311212990, 0.4212463955, 1.2738041E-2, 
                   6.2395678E-2, 6.2395792E-2, 7.2423939E-2, 7.4307981E-2 } };
      principalScatterSAB = sab;
      temps  = { 283.6, 550.0, 600.0, 650.0, 800.0 };
      dwpix  = { 1.656474, 7.2675027, 8.6193808, 12.3081033, 18.58051546 };
      dwp1   = { 0.0, 0.0, 0.0, 0.0, 0.0 };
      tempf  = { 1192.460476, 1278.849394, 1301.020228, 1325.615088, 1412.887566 };
      tempf1 = { 0.0, 0.0, 0.0, 0.0, 0.0 };


      double awr = 0.9991673, spr = 20.43608, aws = 15.85751, sps = 3.7939;
      unsigned int numSecondaryScatterers = 1, secondaryScatterType = 0,
                        numPrincipalAtoms = 2,    numSecondaryAtoms = 1;
      int iel = 0;
      double translationalWeight = 6.9e-3;
      std::vector<double> bragg {};
      int numEdges = 0;
      int za = 1001;
      int ilog = 0, isym = 0, lat = 1; 
    
      std::vector<double> awrVec {awr, aws};
      std::vector<unsigned int> numAtomsVec {numPrincipalAtoms,numSecondaryAtoms};
  

      njoy::ENDFtk::file::Type<7> myMF7 = endout(sab,za,awrVec,spr,sps,temps,
      numSecondaryScatterers,secondaryScatterType,principalScatterSAB,alphas,betas,
      dwpix,dwp1,iel,translationalWeight,bragg,numEdges,tempf,tempf1,ilog,isym,lat,
      numAtomsVec);
  
      auto begin = H2O_1.begin(), end = H2O_1.end();
      long lineNumber = 1;
      StructureDivision division( begin, end, lineNumber );
      njoy::ENDFtk::file::Type<7> goodMF7( division, begin, end, lineNumber );
  
      //checkFullMF7( myMF7, goodMF7, betas );
  
    } // WHEN
  } // GIVEN



  GIVEN( "Be-metal example from ENDFB-VIII.0 library" ){
    WHEN( "Coherent Elastic and Inelastic Scattering" ){
      alphas = { 1e-5, 1e-4, 1e-2, 1.0, 5.0, 10.0, 20.0 };
      betas  = { 0.0 , 1e-4, 1e-2, 1.0, 5.0, 10.0, 20.0 };
      sab    = {{ 7.0370229E-7, 7.0372301E-7, 7.0577387E-7, 1.4116065E-6, 
      1.2170307E-11, 1.7831318E-21, 0.0, 7.0372343E-6, 7.0374415E-6, 7.0579502E-6, 
      1.4115607E-5, 1.2169519E-9, 1.9790182E-18, 0.0, 7.0601224E-4, 7.0603296E-4, 
      7.0808382E-4, 1.4065297E-3, 1.2083139E-5, 2.3385738E-11, 4.0799692E-22, 
      0.0699738674, 0.0699756356, 0.0701506867, 0.0963893250, 0.0609049689, 
      1.1474717E-3, 6.2785457E-8, 0.0831645662, 0.0831667740, 0.0833853408, 
      0.1008018396, 0.1602491842, 0.0722414864, 1.7231248E-3, 0.0307479891, 
      0.0307488804, 0.0308371224, 0.0400176130, 0.0859075534, 0.1082226866, 
      0.0293209593, 5.2427266E-3, 5.2428810E-3, 5.2581677E-3, 6.9719489E-3, 
      0.0184325047, 0.0429634563, 0.0762340185}, 
             {4.0533651E-6, 4.0534147E-6, 
      4.0583290E-6, 6.9324270E-6, 8.2963850E-11,1.5871830E-20, 0.0, 4.0538456E-5, 
      4.0538952E-5, 4.0588092E-5, 6.9320726E-5, 8.2952306E-9, 1.9675371E-17, 0.0, 
      4.1053480E-3, 4.1053972E-3, 4.1102672E-3, 6.8930346E-3, 8.1693568E-5, 
      4.3226337E-10,2.9047298E-20,0.3503790866, 0.3503816671, 0.35063713329, 
      0.3521927490, 0.2108203648, 0.0115059177, 4.1423910E-6, 0.1941325368, 
      0.1941348829, 0.1943671448, 0.2156022593, 0.2556926558, 0.1793803211, 
      0.0190373489, 0.1000797079, 0.1000809325, 0.1002021626, 0.1123720386, 
      0.1569720963, 0.1790575095, 0.0948108014, 0.0383713594, 0.0383718292, 
      0.0384183449, 0.0432273443, 0.0653505809, 0.0950446341, 0.126090255 }};
      principalScatterSAB = {{}};
      temps  = { 500.0, 1200.0 };
      dwpix  = { 1.4490209, 7.58479373 };
      dwp1   = { 0.0, 0.0 };
      tempf  = { 586.94826, 1237.45113 };
      tempf1 = { 0.0, 0.0 };
      double awr = 8.93478, spr = 6.153875, aws = 0.0, sps = 0;
      unsigned int numSecondaryScatterers = 0, secondaryScatterType = 0,
                        numPrincipalAtoms = 1,    numSecondaryAtoms = 0;
      int iel = 2;
      double translationalWeight = 0.0;
      std::vector<double> bragg = Be_metal_bragg;
      int numEdges = 254;
      int za = 126;
      int ilog = 0, isym = 0, lat = 1; 
      int ncold = 0;
      int isabt = 0;
    
      std::vector<double> awrVec {awr};
      std::vector<unsigned int> numAtomsVec {numPrincipalAtoms};
  


      nlohmann::json jsonInput = {
        "npr"   : numPrincipalAtoms,
        "nss"   : numSecondaryScatterers,
        "sps"   : numSecondaryAtoms,
        "ncold" : ncold,
        "isabt" : isabt,
        "lat"   : lat,
        "iel"   : iel, 
        "b7"    : secondaryScattererType }


      njoy::ENDFtk::file::Type<7> myMF7 = endout(sab,za,awrVec,spr,sps,temps,
      numSecondaryScatterers,secondaryScatterType,principalScatterSAB,alphas,betas,
      dwpix,dwp1,iel,translationalWeight,bragg,numEdges,tempf,tempf1,ilog,isym,lat,
      numAtomsVec);
  
      auto begin = BeMetal_1.begin(), end = BeMetal_1.end();
      long lineNumber = 1;
      StructureDivision division( begin, end, lineNumber );
      njoy::ENDFtk::file::Type<7> goodMF7( division, begin, end, lineNumber );
  
      checkFullMF7( myMF7, goodMF7, betas );
  
    } // WHEN
  } // GIVEN


  GIVEN( "SiO2 example from ENDFB-VIII.0 library" ){
    WHEN( "Incoherent Elastic and Inelastic" ){
      alphas = { 0.4, 4.0 };
      betas  = { 0.0, 2.0, 4.0, 8.0, 16.0 };
      sab    = {{ 0.202037135, 0.118385279, 0.019020055, 2.1675111E-3, 4.4194542E-6, 
                 0.029581729, 0.062124313, 0.084467155, 0.0746431711, 0.0147378131 }},
      principalScatterSAB = {{ 0.18326071, 5.09481399E-2, 4.83547456E-2, 7.8537767E-4, 
        4.93779416E-7, 9.310869479E-2, 0.118738373, 0.1007432351, 4.598883932E-2, 
        2.755244914E-3}};
      temps  = { 293.6 };
      dwpix  = { 2.063703012 };
      dwp1   = { 1.893758364 };
      tempf1 = { 486.4483635 };
      tempf  = { 508.3412436 };
      double awr = 27.84423, spr = 2.021, aws = 15.862, sps = 7.4975;
      unsigned int numSecondaryScatterers = 1, secondaryScatterType = 0,
                        numPrincipalAtoms = 1,    numSecondaryAtoms = 1;
      int iel = 0;
      double translationalWeight = 0.0;
      std::vector<double> bragg(0);
      int numEdges = 0;
      int za = 147;
      int ilog = 0, isym = 0, lat = 1; 
    
      std::vector<double> awrVec {awr,aws};
      std::vector<unsigned int> numAtomsVec {numPrincipalAtoms, numSecondaryAtoms};
  
      njoy::ENDFtk::file::Type<7> myMF7 = endout(sab,za,awrVec,spr,sps,temps,
      numSecondaryScatterers,secondaryScatterType,principalScatterSAB,alphas,betas,
      dwpix,dwp1,iel,translationalWeight,bragg,numEdges,tempf,tempf1,ilog,isym,lat,
      numAtomsVec);

      auto begin = SiO2_1.begin(), end = SiO2_1.end();
      long lineNumber = 1;
      StructureDivision division( begin, end, lineNumber );
      njoy::ENDFtk::file::Type<7> goodMF7( division, begin, end, lineNumber );

      checkFullMF7( myMF7, goodMF7, betas );

    } // WHEN
  } // GIVEN
} // TEST CASE
*/
