#include "catch.hpp"
#include "endout/endout_util/inelasticOutput.h"
#include "generalTools/testing.h"


TEST_CASE( "processing inelastic scattering data" ){
  GIVEN( "" ){  

    std::vector<double>
      alphas {0.1, 0.2, 0.3},
      betas  {0.0, 0.2, 0.4, 0.6},
      temps  {296.0, 400.0, 1200.0};
    std::vector<std::vector<double>> 
      fullSAB {
      {3.8033571642194675E-2, 1.3097977662965868E-2, 7.2378904308102836E-003, 7.3984644967871882E-003, 7.2833257288873771E-2, 2.5279372139270825E-2, 1.4062459269729585E-2, 1.4368655814228186E-2, 0.10460945860083895, 3.6593469618555065E-2, 2.0491057645860102E-2, 2.0929367457229134E-2},
      {6.8887759608432159E-2, 2.3220815840088310E-2, 1.2566434914915787E-2, 1.2557653144805678E-2, 0.13084592046566734, 4.4606193587369738E-2, 2.4367346375302308E-2, 2.4330193424512068E-2, 0.18641224657295585, 6.4268300314158977E-2, 3.5435625904678573E-2, 3.5354623207321897E-2},
      {0.57443547123231475, 0.18938839472656757,  9.9835746588985447E-2, 9.5084600390500540E-2, 1.0117147679949494, 0.34596257517155315, 0.18769688181160046, 0.17821769135186777,  1.3376224712991214, 0.47421404450597904, 0.26448395355499854, 0.25053326406628923}};
    int isym = 0, ilog = 0, lat = 1;
    std::cout.precision(15);
    //inelasticOutput(alphas,betas,fullSAB,temps,isym,ilog,lat);
    auto out = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
    auto betaVal = std::get<0>(out);
    auto toWrite = std::get<1>(out);

    std::vector<std::vector<double>> correctOut {{0.03803357, 0.07283326, 0.1046095}, {0.06888776, 0.1308459, 0.1864122}, {0.5744355, 1.011715, 1.337622}};
    std::cout << (toWrite[2]|ranges::view::all) << std::endl;
    std::cout << (correctOut[2]|ranges::view::all) << std::endl;
    REQUIRE( ranges::equal(toWrite[0],correctOut[0],equal) );
    REQUIRE( ranges::equal(toWrite[1],correctOut[1],equal) );
    REQUIRE( ranges::equal(toWrite[2],correctOut[2],equal) );

    


    THEN( "" ){
    } // THEN
  } // GIVEN
} // TEST CASE

