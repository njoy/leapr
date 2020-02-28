#include "catch.hpp"
#include "endout/endout_util/inelasticOutput.h"
#include "generalTools/testing.h"


TEST_CASE( "processing inelastic scattering data" ){
  GIVEN( "" ){  
    std::cout.precision(15);
    int lat = 1;

    std::vector<double>
      alphas {0.1, 0.2, 0.3},
      betas  {0.0, 0.2, 0.4, 0.6},
      temps  {296.0, 400.0, 1200.0};
    std::vector<std::vector<double>> fullSAB {
     { 3.80335716E-2, 1.30979776E-2, 7.237890430E-3, 7.39846449E-3, 7.283325728E-2, 
       2.52793721E-2, 1.40624592E-2, 1.436865581E-2, 0.104609458,   3.659346961E-2, 
       2.04910576E-2, 2.092936745E-2 },
     { 6.88877596E-2, 2.32208158E-2, 1.256643491E-2, 1.255765314E-2, 0.1308459204, 
       4.460619358E-2, 2.436734637E-2, 2.433019342E-2, 0.186412246, 6.426830031E-2, 
       3.543562590E-2, 3.535462320E-2 },
     { 0.5744354712, 0.1893883947,  9.98357465889E-2, 9.508460039E-2, 1.011714767, 
       0.3459625751, 0.1876968818, 0.178217691,  1.3376224712, 0.474214044, 
       0.264483953, 0.2505332640 } };
    
    WHEN( "No cold option used" ){
      AND_WHEN( "Symmetric scattering law requested" ){
        int isym = 0;

        AND_WHEN( "No log option is used" ){
          int ilog = 0;
          auto out = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
          auto betaVal = std::get<0>(out);
          auto toWrite = std::get<1>(out);
      
          std::vector<std::vector<double>> correctOut 
          {{0.03803357, 0.07283326, 0.1046095}, 
           {0.06888776, 0.1308459,  0.1864122}, 
           {0.5744355,  1.011715,   1.337622}};
          REQUIRE( 0.0 == Approx(betaVal).epsilon(1e-6) );
          REQUIRE( ranges::equal(toWrite[0],correctOut[0],equal) );
          REQUIRE( ranges::equal(toWrite[1],correctOut[1],equal) );
          REQUIRE( ranges::equal(toWrite[2],correctOut[2],equal) );
        } // AND WHEN
        AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
          int ilog = 1;
          auto out = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
          auto betaVal = std::get<0>(out);
          auto toWrite = std::get<1>(out);
      
          std::vector<std::vector<double>> correctOut 
           {{-3.269286 -2.619583 -2.257521},
            {-2.675277, -2.033735, -1.679795},
            {-0.5543675,  0.01164668,  0.2908938}};

          REQUIRE( 0.0 == Approx(betaVal).epsilon(1e-6) );
          //REQUIRE( ranges::equal(toWrite[0],correctOut[0],equal) );
          //REQUIRE( ranges::equal(toWrite[1],correctOut[1],equal) );
          //REQUIRE( ranges::equal(toWrite[2],correctOut[2],equal) );

        } // AND WHEN
      } // AND WHEN
      AND_WHEN( "Non-symmetric scattering law requested" ){
        int isym = 2;

        AND_WHEN( "No log option is used" ){
          int ilog = 0;
          auto out = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
          auto betaVal = std::get<0>(out);
          auto toWrite = std::get<1>(out);
    
          std::vector<std::vector<double>> correctOut 
          {{0.03803357, 0.07283326, 0.1046095}, 
           {0.06888776, 0.1308459,  0.1864122}, 
           {0.5744355,  1.011715,   1.337622}};
          REQUIRE( 0.0 == Approx(betaVal).epsilon(1e-6) );
          REQUIRE( ranges::equal(toWrite[0],correctOut[0],equal) );
          REQUIRE( ranges::equal(toWrite[1],correctOut[1],equal) );
          REQUIRE( ranges::equal(toWrite[2],correctOut[2],equal) );
        } // AND WHEN

        AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        } // AND WHEN

      } // AND WHEN

    } // AND WHEN
    THEN( "" ){
    } // THEN
  } // GIVEN
} // TEST CASE

        //inelasticOutput(alphas,betas,fullSAB,temps,isym,ilog,lat);
