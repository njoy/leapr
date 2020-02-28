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
    std::vector<std::vector<double>> correctOut_beta0 (betas.size());
    std::vector<std::vector<double>> correctOut_beta1 (betas.size());
    
    WHEN( "isym = 0. No cold option used and symmetric S(a,b) requested" ){
      int isym = 0;
      AND_WHEN( "No log option is used" ){
        int ilog = 0;

        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);

        correctOut_beta0 = { {0.03803357, 0.07283326, 0.1046095}, 
                             {0.06888776, 0.1308459,  0.1864122}, 
                             {0.5744355,  1.011715,   1.337622 } };
        correctOut_beta1 = { {0.01186118, 0.02289232, 0.03313806}, 
                             {0.02157749, 0.04144943, 0.05972006}, 
                             {0.1848110,  0.3376009,  0.4627526}};
        auto toWrite0 = std::get<1>(out0),
             toWrite1 = std::get<1>(out1);

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );

        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN
      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0);
        auto toWrite1 = std::get<1>(out1);
    
        correctOut_beta0 = {{-3.269286,  -2.619583, -2.257521},
                            {-2.675277,  -2.033735, -1.679795}, 
                            {-0.5543675, 0.01164668, 0.2908938}};

        correctOut_beta1 = {{-4.434485, -3.776954, -3.407073 }, 
                            {-3.836105, -3.183281, -2.818087 }, 
                            {-1.688422, -1.085891, -0.7705627 }};

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN
    } // AND WHEN

    WHEN( "isym = 2. No cold option used and non-symmetric S(a,b) requested" ){
      int isym = 2;

      AND_WHEN( "No log option is used" ){
        int ilog = 0;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0),
             toWrite1 = std::get<1>(out1);
  
        correctOut_beta0 =
        {{0.03803357, 0.07283326, 0.1046095}, 
         {0.06888776, 0.1308459,  0.1864122}, 
         {0.5744355,  1.011715,   1.337622}};
        correctOut_beta1 =
        {{0.01309798, 0.02527937, 0.03659347},
         {0.02322082, 0.04460619, 0.06426830},
         {0.1893884, 0.3459626, 0.4742140}};

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN

      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0);
        auto toWrite1 = std::get<1>(out1);
    
        correctOut_beta0 = {{-3.269286, -2.619583, -2.257521},
                            {-2.675277, -2.033735, -1.679795}, 
                            {-0.5543675, 0.01164668, 0.2908938}};
        correctOut_beta1 = {{-4.335297, -3.677767, -3.307885},
{-3.762706, -3.109883, -2.744689},
{-1.663955, -1.061425, -0.7460965}};

std::cout << (toWrite1[0]|ranges::view::all) << std::endl;
        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          //REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }

      } // AND WHEN

    } // AND WHEN
    THEN( "" ){
    } // THEN
  } // GIVEN
} // TEST CASE

        //inelasticOutput(alphas,betas,fullSAB,temps,isym,ilog,lat);
